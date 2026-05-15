import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from scipy.stats import pearsonr

class BaseEstimator:
    def estimate(self, X, Y, T, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")


class HAIPWEstimator(BaseEstimator):
    def __init__(self, alpha_ridge=0.1, n_folds=5):
        self.alpha_ridge = alpha_ridge
        self.n_folds = n_folds

    def compute_aipw(self, X, Y, T):
        ## CHANGES FROM ORIGINAL REPO: 
        # - no feature selection (not needed for small number of features) 

        kf = KFold(n_splits=self.n_folds, shuffle=True, random_state=42)
        psi_aipw = np.zeros(X.shape[0])
        pi = T.mean() # treatment probability known in the data

        for train_idx, test_idx in kf.split(X):
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            Y_train, Y_test = Y.iloc[train_idx], Y.iloc[test_idx]
            T_train, T_test = T.iloc[train_idx], T.iloc[test_idx]

            mu1_model = Ridge(alpha=self.alpha_ridge)
            mu0_model = Ridge(alpha=self.alpha_ridge)

            mu1_model.fit(X_train[T_train == 1], Y_train[T_train == 1])
            mu0_model.fit(X_train[T_train == 0], Y_train[T_train == 0])

            Q1 = mu1_model.predict(X_test)
            Q0 = mu0_model.predict(X_test)

            psi_aipw[test_idx] = Q1 - Q0 + (T_test * (Y_test - Q1)) / pi - ((1 - T_test) * (Y_test - Q0)) / (1 - pi)

        return psi_aipw

    def estimate(self, X, Y, T, y1_preds, y0_preds, psi_aipw=None):
        # Compute AIPW if not provided
        if psi_aipw is None:
            psi_aipw = self.compute_aipw(X, Y, T)

        pi = T.mean()

        psi_model_estimates = np.array([
            (T * (Y - y1_hat)) / pi - ((1 - T) * (Y - y0_hat)) / (1 - pi) + (y1_hat - y0_hat)
            for y1_hat, y0_hat in zip(y1_preds, y0_preds)
        ])

        Sigma = np.cov(np.vstack((psi_aipw, psi_model_estimates)))
        lambda_star = self.compute_lambda(Sigma)

        haipw_if = (np.vstack((psi_aipw, psi_model_estimates)) * lambda_star[:, np.newaxis]).sum(axis=0)

        return np.mean(haipw_if), lambda_star.T @ Sigma @ lambda_star, lambda_star # 2nd term: matrix product (gives variance)

    def compute_lambda(self, Sigma):
        try:
            Sigma_inv = np.linalg.inv(Sigma)
        except np.linalg.LinAlgError:
            Sigma += 1e-6 * np.eye(Sigma.shape[0])
        Sigma_inv = np.linalg.inv(Sigma)
        ones = np.ones(Sigma.shape[0])
        lambda_star = Sigma_inv @ ones / (ones.T @ Sigma_inv @ ones)
        return lambda_star


class AIPWEstimator(BaseEstimator):
    def __init__(self, alpha_ridge=0.1, n_folds=5):
        self.alpha_ridge = alpha_ridge
        self.n_folds = n_folds

    def estimate(self, X, Y, T, *args, **kwargs):
        kf = KFold(n_splits=self.n_folds, shuffle=True, random_state=42)
        psi_rct = np.zeros(X.shape[0])
        pi = T.mean() # treatment probability known in the data

        for train_idx, test_idx in kf.split(X):
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            Y_train, Y_test = Y.iloc[train_idx], Y.iloc[test_idx]
            T_train, T_test = T.iloc[train_idx], T.iloc[test_idx]

            mu1_model = Ridge(alpha=self.alpha_ridge)
            mu0_model = Ridge(alpha=self.alpha_ridge)

            mu1_model.fit(X_train[T_train == 1], Y_train[T_train == 1])
            mu0_model.fit(X_train[T_train == 0], Y_train[T_train == 0])

            Q1 = mu1_model.predict(X_test)
            Q0 = mu0_model.predict(X_test)

            psi_rct[test_idx] = Q1 - Q0 + (T_test * (Y_test - Q1)) / pi - ((1 - T_test) * (Y_test - Q0)) / (1 - pi)

        return np.mean(psi_rct), np.var(psi_rct, ddof=1)


class PPIEstimator(BaseEstimator):
    def estimate(self, X, Y, T, f_model):
        sigma_t = np.std(Y[T == 1], ddof=1)
        sigma_c = np.std(Y[T == 0], ddof=1)
        sigma_f = np.std(f_model, ddof=1)

        pi_t = T.mean()
        pi_c = 1 - pi_t

        rho_c = pearsonr(f_model[T == 0], Y[T == 0])[0] if np.sum(T == 0) > 0 else 0.0
        rho_t = pearsonr(f_model[T == 1], Y[T == 1])[0] if np.sum(T == 1) > 0 else 0.0

        optimal_lambda = (pi_c * sigma_t * rho_t + pi_t * sigma_c * rho_c) / (sigma_f + 1e-3)
        ppi_est = np.mean(Y[T == 1] - optimal_lambda * f_model[T == 1]) - np.mean(Y[T == 0] - optimal_lambda * f_model[T == 0])

        sigma_t_sq = np.var(Y[T == 1], ddof=1)
        sigma_c_sq = np.var(Y[T == 0], ddof=1)
        sigma_f_sq = np.var(f_model, ddof=1)

        ppi_var = sigma_t_sq / pi_t + sigma_c_sq / pi_c - optimal_lambda**2 * sigma_f_sq * (1 / pi_c + 1 / pi_t)

        return ppi_est, ppi_var

class PPIEstimatorCrossFit(BaseEstimator):
    def estimate(self, X1, Y1, T1, X2, Y2, T2, f_model1, f_model2):

        sigma_t_1 = np.std(Y1[T1 == 1], ddof=1)
        sigma_c_1 = np.std(Y1[T1 == 0], ddof=1)
        sigma_f_1 = np.std(f_model1, ddof=1)
        pi_t_1 = T1.mean()
        pi_c_1 = 1 - pi_t_1
        rho_c_1 = pearsonr(f_model1[T1 == 0], Y1[T1 == 0])[0] if np.sum(T1 == 0) > 0 else 0.0
        rho_t_1 = pearsonr(f_model1[T1 == 1], Y1[T1 == 1])[0] if np.sum(T1 == 1) > 0 else 0.0

        sigma_t_2 = np.std(Y2[T2 == 1], ddof=1)
        sigma_c_2 = np.std(Y2[T2 == 0], ddof=1)
        sigma_f_2 = np.std(f_model2, ddof=1)
        pi_t_2 = T2.mean()
        pi_c_2 = 1 - pi_t_2
        rho_c_2 = pearsonr(f_model2[T2 == 0], Y2[T2 == 0])[0] if np.sum(T2 == 0) > 0 else 0.0
        rho_t_2 = pearsonr(f_model2[T2 == 1], Y2[T2 == 1])[0] if np.sum(T2 == 1) > 0 else 0.0

        optimal_lambda_1 = (pi_c_1 * sigma_t_1 * rho_t_1 + pi_t_1 * sigma_c_1 * rho_c_1) / (sigma_f_1 + 1e-3)
        optimal_lambda_2 = (pi_c_2 * sigma_t_2 * rho_t_2 + pi_t_2 * sigma_c_2 * rho_c_2) / (sigma_f_2 + 1e-3)
        ppi_est_1 = np.mean(Y1[T1 == 1] - optimal_lambda_2 * f_model1[T1 == 1]) - np.mean(Y1[T1 == 0] - optimal_lambda_2 * f_model1[T1 == 0])
        ppi_est_2 = np.mean(Y2[T2 == 1] - optimal_lambda_1 * f_model2[T2 == 1]) - np.mean(Y2[T2 == 0] - optimal_lambda_1 * f_model2[T2 == 0])

        ppi_est = 0.5 * (ppi_est_1 + ppi_est_2)
        return ppi_est

class DifferenceInMeansEstimator(BaseEstimator):
    def estimate(self, X, Y, T, *args, **kwargs):
        pi = T.mean()
        dm_i = Y * T / pi - Y * (1 - T) / (1 - pi)
        dm_est = dm_i.mean()
        dm_var = np.var(Y[T == 1], ddof=1) / pi + np.var(Y[T == 0], ddof=1) / (1 - pi)

        return dm_est, dm_var
