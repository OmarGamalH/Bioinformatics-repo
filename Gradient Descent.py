import numpy as np
def y_pred(X , W , B):
    return np.dot(X , W) + B


def loss_func(X , Y , W , B , y_pred):
    """
    X : MATRIX
    Y : MATRIX
    W : VECTOR
    B : SCALAR
    y_pred : function
    
    
    """
    m , n = X.shape
    cost = 0
    for i in range(m):
        cost = cost + (y_pred(X[i] , W , B) - Y[i]) ** 2
    
    cost = cost /(2 * m)

    return cost


def compute_gradient(X , Y , W , B , y_pred):
    """
    X : MATRIX
    Y : MATRIX
    W : VECTOR
    B : SCALAR
    y_pred : function
    learning_rate : scalar
    
    """
    m , n = X.shape
    
    dl_dw = np.zeros(n)
    error = 0
    for i in range(m):
        err_i = y_pred(X[i] , W , B) - Y[i]
        for j in range(n): # j = 1
            dl_dw[j] = dl_dw[j] + err_i[0] * X[i , j] # dl_dw[1] = dl_dw[1] + err_i * X[i , 1]
        error = error + err_i

    dl_db = (1/m) * error
    dl_dw = dl_dw * (1/m)

    return dl_db[0] , dl_dw
    


def gradient_descent(X , Y , W , B , iterations ,  y_pred , compute_gradient , learning_rate):
    W_init = W.copy()
    iters = []
    losses = []
    print(f"X:{X.shape} , Y:{Y.shape} , W:{W.shape} , B:{B}")
    for i in range(iterations):
        dl_db , dl_dw = compute_gradient(X , Y , W_init , B , y_pred)
        print(f"W_init:{W_init} , B:{B}")
        losses.append(loss_func(X , Y , W_init , B , y_pred))
        iters.append(i + 1)
        W_init -= learning_rate * dl_dw
        B -= learning_rate * dl_db
        
    
    return W_init , B , losses , iters 


# W, W_B , losses , iters = gradient_descent(X_TRAIN , Y_TRAIN , W_init , B , 100 , y_pred , compute_gradient , learning_rate=1)