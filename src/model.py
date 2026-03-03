import pandas as pd
import numpy as np
from tabpfn import TabPFNRegressor
import argparse
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

FEATURE_COL_START = 7
FEATURE_COL_END = 191
TARGET_COL = 'Coupling_Value'
OUTPUT_COLS = ['GPCR', 'Gprotein_name', 'Coupling_Value']
MIN_COUPLING_VALUE = 0.01

def load_and_preprocess_training(file_path):
    """Loads training data and preprocesses the target variable."""
    known = pd.read_feather(file_path)

    known[TARGET_COL] = known[TARGET_COL].apply(lambda x: max(MIN_COUPLING_VALUE, x))
    known[TARGET_COL] = np.log10(known[TARGET_COL])

    X = known.iloc[:, FEATURE_COL_START:FEATURE_COL_END]
    y = known[TARGET_COL]

    feature_names = X.columns.tolist()

    return X, y, feature_names

def train_model(X_train, y_train, device):
    """Initializes and trains the TabPFN model."""
    model = TabPFNRegressor(device=device)
    model.fit(X_train, y_train)
    return model

def predict_and_postprocess(model, input_file, feature_names):
    """Loads data, predicts, post-processes, and returns results."""
    data = pd.read_csv(input_file, sep="\t")
    X_data = data[feature_names]
    predictions_log = model.predict(X_data)
    predictions = np.power(10, predictions_log)
    data[TARGET_COL] = predictions
    return data

def save_results(data, output_file):
    """Saves the specified columns to the output file."""
    print(f"Saving results to: {output_file}")
    output_df = data[OUTPUT_COLS]
    output_df.to_csv(output_file, sep="\t", index=False)
    return

def main():
    """Main function to orchestrate loading, training, prediction, and saving."""
    parser = argparse.ArgumentParser(description="Train a TabPFN model and predict coupling values.")

    parser.add_argument("--input_file",
                        help="Path to the input TSV file for prediction.")
    parser.add_argument("--output_file",
                        help="Path to save the output TSV file with predictions.")
    parser.add_argument("--training-data",
                        default="data/known_pairs.feather",
                        help="Path to the training data Feather file")
    parser.add_argument("--device", "-d",
                        default="cuda:0",
                        help="Device to use for TabPFN")

    args = parser.parse_args()

    X_train, y_train, feature_names = load_and_preprocess_training(args.training_data)
    model = train_model(X_train, y_train, args.device)
    results_data = predict_and_postprocess(model, args.input_file, feature_names)
    save_results(results_data, args.output_file)

if __name__ == "__main__":
    main()