import numpy as np 
import pandas as pd
from tdc.multi_pred import DrugSyn
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split, GridSearchCV
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_absolute_error as mae, mean_squared_error as mse
from scipy import stats

def load_data(data_path, data_file):
	df = pd.read_csv(data_path + data_file + ".csv")
	return df

def preprocess_data(df):
	df = df.dropna(subset=['AAC']) #Change IC50 for AAC, depending on the desired metric
	X = df.drop(columns=['cell_line', 'AAC', 'IC50']) #Change IC50 for AAC
	y = df['AAC'] #Change IC50 for AAC
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
	data_path = "./data/"
	data_file = "merged_GDSC" #change to merged_NCI for NCI data
	
	df = load_data(data_path, data_file)
	print(df.shape)
	print(df.columns)
	
	X, y = preprocess_data(df)
	X_train, X_test, y_train, y_test = split_data(X, y)
	
	clf = train_model(X_train, y_train)
	
	corr, abs_error, square_error, preds = evaluate_model(clf, X_test, y_test)
	print(f"\n\tThe correlation value between the predicted and true synergy values is {np.round(corr, 3)}")
	print(f"\n\tThe mean absolute error is {np.round(abs_error, 3)}")
	print(f"\n\tThe mean square error is {np.round(square_error, 3)}")
	
	visualize_results(y_test, preds)

if __name__ == "__main__":
	main()
