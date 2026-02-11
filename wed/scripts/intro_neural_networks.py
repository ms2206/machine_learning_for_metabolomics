# import libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tensorflow.keras.datasets import mnist
from tensorflow.keras.utils import to_categorical
from tensorflow import keras
from tensorflow.keras import layers

# load the dataset
(x_train, y_train), (x_test, y_test) = mnist.load_data()

# Flatten the images
x_train = x_train.reshape(-1, 784)
x_test = x_test.reshape(-1, 784)

# standardize the pixel values by dividing by 255
x_train = x_train.astype('float32') / 255.0
x_test = x_test.astype('float32') / 255.0

# convert the vecotor containing the train set into a matrix with 10 columns (one for each class)
y_train_one_hot = to_categorical(y_train, num_classes=10)
y_test_one_hot = to_categorical(y_test, num_classes=10)

# change the column names for the mnist matrix to be between 0 and 9
column_names = [str(i) for i in range(10)]
y_train_one_hot_df = pd.DataFrame(y_train_one_hot, columns=column_names)

# print(y_train_one_hot_df.head())

# define the model architecture
model = keras.Sequential([
    keras.Input(shape=(784,)),
    layers.Dense(256, activation='relu'),
    layers.Dense(128, activation='relu'),
    layers.Dense(64, activation='relu'),
    layers.Dense(10, activation='softmax')
])

# print the model summary
# model.summary()

# compile the model
model.compile(
    loss='categorical_crossentropy',
    optimizer=keras.optimizers.RMSprop(),
    metrics=['accuracy']
)

# train the model
history_1 = model.fit(
    x_train, 
    y_train_one_hot,
    epochs=30,
    batch_size=128,
    validation_split=0.2,
    verbose=0)

# print(history_1.history.keys())

