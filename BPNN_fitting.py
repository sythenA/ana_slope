
import numpy as np
from keras.models import Sequential
from keras.datasets import mnist
from keras.layers import Dense, Activation
from matplotlib import pyplot as plt


def BPNN_model(Input, Output):
    model = Sequential()
    model.add(Dense(units=200, input_dim=2, kernel_initializer='normal',
                    activation='relu'))
    model.add(Dense(units=20, kernel_initializer='normal',
                    activation='relu'))
    model.add(Dense(units=2, kernel_initializer='normal',
                    activation='relu'))
    model.compile(loss='mean_squared_error', optimizer='sgd',
                  metrics=['accuracy'])
    train_history = model.fit(x=Input, y=Output,
                              epochs=150, verbose=1)
    return model
