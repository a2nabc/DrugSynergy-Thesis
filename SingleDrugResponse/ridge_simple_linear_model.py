import numpy as np 
import pandas as pd
from tdc.multi_pred import DrugSyn
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split, GridSearchCV
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from sklearn.metrics import mean_absolute_error as mae, mean_squared_error as mse
from scipy import stats
import os

def load_data(data_path, data_file):
	df = pd.read_csv(data_path + data_file + ".csv")
	return df

def preprocess_data(df, metric):
	if metric == "AAC":
		df = df.drop(columns=['IC50'])
	else:
		df = df.drop(columns=['AAC'])
	df = df.dropna(subset=[metric])
	X = df.drop(columns=['cell_line', metric])
	y = df[metric]
	return X, y

def split_data(X, y, train_size=0.8, test_size=0.2, seed=1234):
	return train_test_split(X, y, train_size=train_size, test_size=test_size, random_state=seed)

def train_model(X_train, y_train, num_folds=5):
	model = Ridge()
	parameters = {'alpha': [0.05, 0.001, 0.1, 0.5, 1, 2, 10]}
	clf = GridSearchCV(model, parameters, cv=num_folds)
	clf.fit(X_train, y_train)
	return clf

def evaluate_model(clf, X_test, y_test):
	preds = clf.predict(X_test)
	corr, _ = stats.pearsonr(y_test, preds)
	abs_error = mae(y_test, preds)
	square_error = mse(y_test, preds)
	return corr, abs_error, square_error, preds

def visualize_results(y_test, preds):
	res_df = pd.DataFrame({'True Value': y_test, 'Predicted Value': preds})
	sns.scatterplot(data=res_df, x="True Value", y="Predicted Value")
	plt.show()


def main():
	data_path = "./data/allGeneExpressions/"
	data_files = [f.split(".csv")[0] for f in os.listdir(data_path) if f.endswith(".csv")]

	metrics = ["AAC", "IC50"]
	for metric in metrics:
		results = []
		print(f"Running for metric: {metric}")

		for data_file in data_files:
			df = load_data(data_path, data_file)
			print(f"Processing {data_file}: {df.shape}")
			print(df.columns)

			X, y = preprocess_data(df, metric)

			corrs, abs_errors, square_errors, all_preds, all_true = [], [], [], [], []
			for seed in range(10):
				seed = np.random.randint(0, 10000)
				X_train, X_test, y_train, y_test = split_data(X, y, seed=seed)
				clf = train_model(X_train, y_train)
				corr, abs_error, square_error, preds = evaluate_model(clf, X_test, y_test)
				corrs.append(corr)
				abs_errors.append(abs_error)
				square_errors.append(square_error)
				print(f"File {data_file}, Seed {seed}: Correlation = {np.round(corr, 3)}, MAE = {np.round(abs_error, 3)}, MSE = {np.round(square_error, 3)}")
				all_preds.extend(preds)
				all_true.extend(y_test)

			avg_corr = np.round(np.mean(corrs), 3)
			avg_mae = np.round(np.mean(abs_errors), 3)
			avg_mse = np.round(np.mean(square_errors), 3)
			results.append([data_file, avg_corr, avg_mae, avg_mse])

			# Visualize and save the results
			visualize_results(all_true, all_preds)
			plt.savefig(f"{data_file}_{metric}_results.png")
			plt.close()

		results_df = pd.DataFrame(results, columns=["File", "Average Correlation", "Average MAE", "Average MSE"])
		print(results_df)
		results_df.to_csv(f"results_summary_{metric}.csv", index=False)

if __name__ == "__main__":
	main()

