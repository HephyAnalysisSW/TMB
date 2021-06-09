from keras.models import load_model

keras_model_loaded = load_model("TTZ_TWZ_WZ_Keras_Model.h5")
pred = keras_model_loaded.predict(X_test[:3])
print(pred)


