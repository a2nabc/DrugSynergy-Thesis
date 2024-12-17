"""
Example code to predict drug synergy using simple linear models.

This code is built on data available from the Therapeutic Data Commons (tdcommons)
and uses the OncoPolyPharmacology screen. tdcommons has another synergy data set available
(DrugComb)but let's ignore that for now. 

Some notes:
 - OncoPolyPharmacology only has the Loewe Additivity synergy score
 - OncoPolyPharmacology is a large dataset (~1.6 GB)


In brief what the code below does is select two drugs and focus on their synergy
values across a number of cell lines. This is very similar to the task you're doing
with GDSC except you have one drug. This should give you some example code to get 
started on running a simpler model. I use a Lasso model to give you a chance to look at
regularization. The wikipedia is decent: https://en.wikipedia.org/wiki/Lasso_(statistics)


The datasets are not the one Ben wants you to use but it's good practice for you
and also can be used as experiments in your thesis. 



Written by: James Bannon




"""

# required libraries

import numpy as np 
import pandas as pd
from tdc.multi_pred import DrugSyn # see https://tdcommons.ai for install instructions
from sklearn.linear_model import LinearRegression, ElasticNet, Lasso, Ridge
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_squared_error as mse
from scipy import stats

data_path = "./data/"
data_NCI = "merged_NCI"
data_GDSC = "merged_GDSC"

# this will automatically download a drug synergy dataset when you first call it.
# subsequent calls won't need a download
# docs: https://tdc.readthedocs.io/en/main/tdc.multi_pred.html

df = pd.read_csv(data_path+data_NCI+".csv")       #Change NCI for GDSC
print(df.shape)
print(df.columns)

# Drop rows with missing IC50 values
df = df.dropna(subset=['IC50']) #change this to compare IC50 vs AAC

# take the cell line gene expression values and make them a matrix
X = df.drop(columns=['cell_line', 'AAC', 'IC50'])
y = df['IC50'] #change this to compare IC50 vs AAC

#print(X.shape)
#print(y.shape)

train_size, test_size = 0.8, 0.2
seed = 1234

X_train, X_test, y_train, y_test = train_test_split(
	X,
	y,
	train_size=train_size,
	test_size=test_size,
	random_state = seed)
# train_test_split(*arrays, test_size=None, train_size=None, random_state=None, shuffle=True, stratify=None)[source]


num_folds = 3
#specify the model
model = Lasso()
#specify the grid of parameters you want to try

parameters = {'alpha':[0.05,0.001,0.1,0.5,1,2,10]}

# see docs: https://scikit-learn.org/1.5/modules/generated/sklearn.model_selection.GridSearchCV.html
clf = GridSearchCV(model, parameters,cv =num_folds) 

clf.fit(X_train,y_train)
preds = clf.predict(X_test) # this will automatically use the best performing hyperparameters.

## evaluating the model:
corr,_= stats.pearsonr(y_test, preds)
abs_error = mae(y_test,preds)
square_error = mse(y_test,preds)
print(f"\n\tThe correlation value between the predicted and true synergy values is {np.round(corr,3)}")
print(f"\n\tThe mean absolute error is {np.round(abs_error,3)}")
print(f"\n\tThe mean square error is {np.round(square_error,3)}")


# Visualize the predicted versus true. Good prediction should look like a straight line.
res_df = pd.DataFrame({'True Value': y_test,'Predicted Value':preds})
sns.scatterplot(data=res_df, x="True Value", y="Predicted Value")
plt.show()


