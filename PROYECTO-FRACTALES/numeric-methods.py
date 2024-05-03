## IMPLEMENTACION DE LOS METODOS

### Newton clasico

def Newton_classicroot(func, dfunc, x0, tol, maxite):
    """
    Implementación del método de Newton para encontrar las raíces de una función en números complejos.
    
    Parámetros:
    - func: la función dada.
    - dfunc: las derivadas de la función f.
    - x0: valor inicial para comenzar la iteración.
    - tol: tolerancia para la convergencia.
    - maxite: número máximo de iteraciones permitidas.
    
    Retorna:
    - Las iteraciones del metodo usando x0 y la raiz encontrada x
    """
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root - x)
                return ite, x
            x = x - func(x)/dfunc(x)[0]
        except Exception:
            ite = maxite # Indicamos que no hay convergencia a partir del x0
            break
    return ite, x

### Newton raices multiples

def Newton_multipleroot(func, dfunc, x0, tol, maxite):
    """
    Implementación del método de Newton para encontrar las raíces complejas multiples.
    
    Parámetros:
    - func: la función dada.
    - dfunc: las derivadas de la función f.
    - x0: valor inicial para comenzar la iteración.
    - tol: tolerancia para la convergencia.
    - maxite: número máximo de iteraciones permitidas.
    
    Retorna:
    - Las iteraciones del metodo usando x0 y la raiz encontrada x
    """    
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root - x)
                return ite, x
            x = x - (func(x)*dfunc(x)[0])/((dfunc(x)[0])**2 - func(x)*dfunc(x)[1])
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

### Aceleracion convexa de Wittaker

def Convex_acc_Wh(func,dfunc,x0,tol,maxite):
    """
    Implementación del método de Newton para encontrar las raíces complejas multiples.
    
    Parámetros:
    - func: la función dada.
    - dfunc: las derivadas de la función f.
    - x0: valor inicial para comenzar la iteración.
    - tol: tolerancia para la convergencia.
    - maxite: número máximo de iteraciones permitidas.
    
    Retorna:
    - Las iteraciones del metodo usando x0 y la raiz encontrada x
    """    
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            d2fx = dfunc(x)[1]
            L = (fx * d2fx) / (dfx)**2
            x = x - (fx/(2*dfx))*(2-L)
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

### Aceleracion doble de Wittaker

def Double_convex_acc_Wh(func,dfunc,x0,tol,maxite):
    """
    Implementación del método de Newton para encontrar las raíces complejas multiples.
    
    Parámetros:
    - func: la función dada.
    - dfunc: las derivadas de la función f.
    - x0: valor inicial para comenzar la iteración.
    - tol: tolerancia para la convergencia.
    - maxite: número máximo de iteraciones permitidas.
    
    Retorna:
    - Las iteraciones del metodo usando x0 y la raiz encontrada x
    """    
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            d2fx = dfunc(x)[1]
            L = (fx * d2fx) / (dfx)**2
            x = x - (fx/(4*dfx))*(2-L + (4+2*L)/(2-L*(2-L)))
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Halley(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            d2fx = dfunc(x)[1]
            L = (fx * d2fx) / (dfx)**2
            x = x - (fx/dfx) * (2/(2 - L))
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Chebyshev(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            d2fx = dfunc(x)[1]
            L = (fx * d2fx) / (dfx)**2
            x = x - (fx/dfx) * (1 + L/2)
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def ConvexN_sHalley(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            d2fx = dfunc(x)[1]
            L = (fx * d2fx) / (dfx)**2
            x = x - (fx/(2*dfx)) * (2 - L)/(1 - L)
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Stir(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            x = x - func(x)/dfunc(x-func(x))[0]
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Steffensen(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            gx = (func(x + func(x)) - func(x))/func(x)
            x = x - fx/gx
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Midpoint_Mth(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            x = x - fx/dfunc(x - fx/(2*dfx))[0]
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Traub_Ostro(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            ux = fx/dfx
            x = x - ux*((func(x-ux)-fx)/(2*func(x-ux)-fx))
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Jarrat(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            ux = fx/dfx
            x = x - 0.5*ux + fx/(dfx - 3*dfunc(x-(2/3)*ux)[0])
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Inverse_Jarrat(func,dfunc,x0,tol,maxite):
    x = x0  
    for ite in range(maxite+1):
        try:
            if abs(func(x)) < tol: # abs(func(x)) abs(root -x)
                return ite, x
            fx = func(x)
            dfx = dfunc(x)[0]
            ux = fx/dfx
            hx = (dfunc(x-(2/3)*ux)[0] - dfx)/dfx  
            x = x - ux + 0.75*ux*hx*(1 - 1.5*hx)
        except Exception:
                ite = maxite # Indicamos que no hay convergencia a partir del x0
                break
    return ite, x

def Numeric_Method(num,func,dfunc,x0,tol,maxite):
    Methods = [Newton_classicroot,Newton_multipleroot,Convex_acc_Wh,Double_convex_acc_Wh,Halley,Chebyshev,ConvexN_sHalley,Stir,Steffensen,Midpoint_Mth,Traub_Ostro,Jarrat,Inverse_Jarrat]
    return Methods[num](func,dfunc,x0,tol,maxite)