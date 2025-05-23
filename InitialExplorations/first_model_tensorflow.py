import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.optimizers import Adam

# Load dataset with synergy and ic50 data
file_path = 'Data/orcestra/synergy_and_ic50_data.csv'
data = pd.read_csv(file_path)

# Encode categorical variables
encoder = LabelEncoder()
data['drugid1_encoded'] = encoder.fit_transform(data['drugid1'])
data['drugid2_encoded'] = encoder.fit_transform(data['drugid2'])

# Features and target
X = data[['drugid1_encoded', 'drugid2_encoded', 'ic50_drug1', 'ic50_drug2']].values
y = data['synergy_score'].values

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Choose between quartile filtering or log scaling
def preprocess_data(X_train, X_test, y_train, y_test, method='quartile'):
    if method == 'quartile':
        # Separate ic50_drug1 and ic50_drug2 into quartiles and eliminate the highest quartile
        ic50_drug1 = X_train[:, 2]
        ic50_drug2 = X_train[:, 3]

        # Calculate the 75th percentile for ic50_drug1 and ic50_drug2
        q75_drug1 = np.percentile(ic50_drug1, 75)
        q75_drug2 = np.percentile(ic50_drug2, 75)

        # Filter out the highest quartile
        mask = (ic50_drug1 <= q75_drug1) & (ic50_drug2 <= q75_drug2)
        X_train = X_train[mask]
        y_train = y_train[mask]

        # Apply the same filtering to the test set
        ic50_drug1_test = X_test[:, 2]
        ic50_drug2_test = X_test[:, 3]
        mask_test = (ic50_drug1_test <= q75_drug1) & (ic50_drug2_test <= q75_drug2)
        X_test = X_test[mask_test]
        y_test = y_test[mask_test]
    elif method == 'log':
        # Apply log scaling to ic50 values
        X_train[:, 2:] = np.log1p(X_train[:, 2:])
        X_test[:, 2:] = np.log1p(X_test[:, 2:])
    else:
        raise ValueError("Method must be either 'quartile' or 'log'")
    
    return X_train, X_test, y_train, y_test

# Preprocess data
X_train, X_test, y_train, y_test = preprocess_data(X_train, X_test, y_train, y_test, method='quartile')

# Normalize IC50 values
scaler = StandardScaler()
X_train[:, 2:] = scaler.fit_transform(X_train[:, 2:])
X_test[:, 2:] = scaler.transform(X_test[:, 2:])

# Build the model
model = Sequential([
    Dense(128, activation='relu', input_shape=(X_train.shape[1],)),
    Dropout(0.3),
    Dense(64, activation='relu'),
    Dense(1)  # Regression output layer
])

# Compile the model
model.compile(optimizer=Adam(learning_rate=0.001), loss='mse', metrics=['mae'])

# Train the model
history = model.fit(X_train, y_train, validation_split=0.2, epochs=50, batch_size=32, verbose=1)

# Evaluate the model
test_loss, test_mae = model.evaluate(X_test, y_test, verbose=1)
print(f"Test Loss: {test_loss}, Test MAE: {test_mae}")

# Predict synergy scores
predicted_synergy = model.predict(X_test)
