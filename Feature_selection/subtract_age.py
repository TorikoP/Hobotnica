import os
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from tqdm import tqdm

def substract_age(file_path):
    data = pd.read_pickle(file_path)
    data_1 = data.copy()
    
    data_1['Condition'] = data_1['Condition'].apply(lambda x: 0 if x == 'HC' else 1)

    for idx, site in tqdm(enumerate(data_1.columns[1:-2])):
        y_train = data_1[site]
        X_train = data_1[['Age', 'Condition']]

        model = LinearRegression()
        model.fit(X_train, y_train)

        age_coef = model.coef_[0]
        intercept = model.intercept_
        data_1[site] = data_1[site] - (age_coef * data_1['Age'] + intercept)
        
        # every 1000th site name for progress tracking
        if (idx + 1) % 1000 == 0:
            print(f'Processed site: {site}')
            
    data_1["Condition"] = data["Condition"]
    
    return data_1

def process_and_save(file_path, output_folder):
    print(f'Starting processing of dataset: {os.path.basename(file_path)}')
    output_path = os.path.join(output_folder, os.path.basename(file_path))

    if os.path.exists(output_path):
        print(f'Skipping already processed file: {output_path}')
        return
    
    #save
    processed_df = substract_age(file_path)
    processed_df.to_pickle(output_path)
    print(f'Processed and saved: {output_path}')

if __name__ == '__main__':
    initial_datasets = '/home/vpal/hobotnica/All_datasets/data_imputed_with_meta_for_regression_pkl'
    substracted_age = '/home/vpal/hobotnica/All_datasets/subtracted_ages_for_all_datasets'

    #
    files = [os.path.join(initial_datasets, filename) for filename in os.listdir(initial_datasets) if filename.endswith('.pickle')]

    for file_path in files:
        process_and_save(file_path, substracted_age)
