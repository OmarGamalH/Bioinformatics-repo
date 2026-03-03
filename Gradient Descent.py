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




class logistics_regression:
    
    def cost_f(self , X , Y , W , B):

        m = X.shape[0]
        cost = 0
        for i in range(m):
            f = np.dot(W , X[i]) + B
            F_X_i = self.sigmoid_f(f)
            loss = Y[i] * np.log(F_X_i) + (1-Y[i]) * np.log(1 - F_X_i)
            cost = cost + loss
        cost = cost / (-1 * m)
        return cost

    def sigmoid_f(self , x):
        y = 1 + np.exp(-1 * x)
        y = 1/y
        return y

    def derivative(self, X , Y , W , B):
        m , n = X.shape

        dj_db = 0
        dj_dw = np.zeros(n)
        for i in range(m):
            f = np.dot(X[i] , W) + B
            F_X_i = self.sigmoid_f(f)
            error = F_X_i - Y[i]
            dj_db += error
            dj_dw += error * X[i]

        dj_db = dj_db / m
        dj_dw = dj_dw / m
        
        return dj_dw , dj_db
    
    def GD_logistic_regression( self, X , Y , W , B , alpha , iterations):

        W_init = W.copy()
        losses = []
        all_iter = []
        for i in range(iterations):
            losses.append(self.cost_f(X , Y , W_init , B))
            all_iter.append(i)
            dj_dw , dj_db = self.derivative(X , Y , W_init , B)
            print(W_init , B)
            W_init -= alpha * dj_dw 
            B -= alpha * dj_db
        
        return W_init , B , losses , all_iter


    def classify(self , row):
        if row >= .5:
            return 1
        return 0
    

    def predict(self , X , W , B):
        sigm = []
        for x in X:
            f = np.dot(W , x) + B 
            sigm.append(self.sigmoid_f(f))
        
        result = [*map(self.classify , sigm)]
        return result