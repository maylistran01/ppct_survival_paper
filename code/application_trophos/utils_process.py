import pandas as pd  
import numpy as np

#################################################################
#### Functions to process dataframes for survival analysis ####
#################################################################


def process_delta_and_filter(df, age_at_death_column="AGE_AT_DEATH", last_known_living_age_column="LAST_KNOWN_LIVING_AGE", covariates=["SEX", "AGE_AT_DIAGNOSIS", "ONSET_SYMPTOMS", "ALSFRS_R_TOTAL", "BMI"]):
    """Prepare event dataframe for survival analysis.

    Args:
        df (pd.DataFrame): Input dataframe with patient data.


    Returns:        
        pd.DataFrame: DataFrame with ID, EVENT, and TIME_OF_EVENT.
    """    

    # Delay to do survival in TIME and not in AGE
    df[f"DELTA_{age_at_death_column}"] = df[age_at_death_column] - df["AGE_AT_BASELINE"]
    df[f"DELTA_{last_known_living_age_column}"] = df[last_known_living_age_column] - df["AGE_AT_BASELINE"]

    # Filter to only keep patients with complete covariates
    df_filtered = df.dropna(subset=covariates)
    return df_filtered

   

def create_event_dataframe(df, delta_age_at_death_col="DELTA_AGE_AT_DEATH", delta_last_known_living_age_col="DELTA_LAST_KNOWN_LIVING_AGE"):
    """Prepare event dataframe for survival analysis.
    Args:
        df (pd.DataFrame): Input dataframe with patient data.
        delta_age_at_death_col (str): Column name for age at death delta.
        delta_last_known_living_age_col (str): Column name for last known living age delta.

    Returns:
        pd.DataFrame: DataFrame with ID, EVENT, and TIME_OF_EVENT.
    """
    # Create event df
    uncensored = df.groupby("ID")[delta_age_at_death_col].first().dropna()
    censored = df[df[delta_age_at_death_col].isna()].groupby("ID")[delta_last_known_living_age_col].first()
    assert(len(censored) + len(uncensored) == df.index.get_level_values("ID").nunique())  # verify that all the patients have events

    uncensored_df = uncensored.reset_index()
    uncensored_df['EVENT'] = True

    censored_df = censored.reset_index()
    censored_df['EVENT'] = False

    uncensored_df = uncensored_df.rename(columns={delta_age_at_death_col: "TIME_OF_EVENT"})
    censored_df = censored_df.rename(columns={delta_last_known_living_age_col: "TIME_OF_EVENT"})
    combined_df = pd.concat([uncensored_df, censored_df], ignore_index=True)

    return combined_df

def simulate_age_at_death(df, intervals):
    df = df.copy()
    df["DELTA_AGE_MIN_AT_DEATH_SIM"] = np.nan

    for i, row in df.iterrows():
        min_age = row["DELTA_AGE_MIN_AT_DEATH"]
        max_age = row["DELTA_AGE_MAX_AT_DEATH"]

        if pd.isna(max_age):  # censored case
            if round(min_age, 6) < 1.5:  # censored during CT
                # find the next interval boundary
                next_interval = next((x for x in intervals if x > min_age), None)
                df.loc[i, "DELTA_LAST_KNOWN_LIVING_AGE_SIM"] = np.random.uniform(min_age, next_interval)
            else:  # censored after CT ended
                df.loc[i, "DELTA_LAST_KNOWN_LIVING_AGE_SIM"] = min_age
        else:  # known death age
            df.loc[i, "DELTA_AGE_AT_DEATH_SIM"] = np.random.uniform(min_age, max_age)

    return df