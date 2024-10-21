import numpy as np


# Xs - HTA - ITA combine
print("============ Xs - HTA - ITA ============")
Covariance = np.array([[2.64**2, 1.01, 0.19], [1.01, 2.3**2, 0.26], [0.19, 0.26, 1.5**2]])
Covariance_inverse = np.linalg.inv(Covariance)  

w1 = np.sum(Covariance_inverse[0, :]) / np.sum(Covariance_inverse)
w2 = np.sum(Covariance_inverse[1, :]) / np.sum(Covariance_inverse)
w3 = np.sum(Covariance_inverse[2, :]) / np.sum(Covariance_inverse)

mu_bar = 1.1*w1+2.2*w2+5.4*w3

chi2 = np.array([[1.1 - mu_bar, 2.2-mu_bar, 5.4-mu_bar]]) @ Covariance_inverse @ np.array([[1.1-mu_bar], [2.2-mu_bar], [5.4-mu_bar]])

print("Covariance matrix:")
print(Covariance)
print("\n")

print("Inverse of covariance matrix:")
print(Covariance_inverse)
print("\n")

print("w1:")
print("%f / %f = %f" % (np.sum(Covariance_inverse[0, :]), np.sum(Covariance_inverse), w1))
print("\n")

print("w2:")
print("%f / %f = %f" % (np.sum(Covariance_inverse[1, :]), np.sum(Covariance_inverse), w2))
print("\n")

print("w3:")
print("%f / %f = %f" % (np.sum(Covariance_inverse[2, :]), np.sum(Covariance_inverse), w3))
print("\n")

print("mu_bar:")
print("%f" % mu_bar)
print("\n")

print("chi2:")
print("%f" % chi2)
print("\n")
print("========================================")


# Xs - HTA combine
print("============ Xs - HTA ==================")
Covariance = np.array([[2.64**2, 1.01], [1.01, 2.3**2]])
Covariance_inverse = np.linalg.inv(Covariance)  

w1 = np.sum(Covariance_inverse[0, :]) / np.sum(Covariance_inverse)
w2 = np.sum(Covariance_inverse[1, :]) / np.sum(Covariance_inverse)

mu_bar = 1.1*w1+2.2*w2

chi2 = np.array([[1.1 - mu_bar, 2.2-mu_bar]]) @ Covariance_inverse @ np.array([[1.1-mu_bar], [2.2-mu_bar]])

print("Covariance matrix:")
print(Covariance)
print("\n")

print("Inverse of covariance matrix:")
print(Covariance_inverse)
print("\n")

print("w1:")
print("%f / %f = %f" % (np.sum(Covariance_inverse[0, :]), np.sum(Covariance_inverse), w1))
print("\n")

print("w2:")
print("%f / %f = %f" % (np.sum(Covariance_inverse[1, :]), np.sum(Covariance_inverse), w2))
print("\n")


print("mu_bar:")
print("%f" % mu_bar)
print("\n")

print("chi2:")
print("%f" % chi2)
print("\n")
print("========================================")
