from sympy.solvers import solve as eq
from sympy import symbols
from numpy import pi, sqrt, arange
import matplotlib.pyplot as plt
import mplcyberpunk
 
 
#Infos
#AÇO AISI 1080
D = 0.9 #m
t = 3.2*10**-3 #m
w = 115 #N/m
l = 4.5 #m
p = [41, 55] #kpa, [min, max]
coef_seg = 3
 
 
#infos tabeladas do livro
sut = 675 #mPa
ai = 0.5 #mm 
 
sy = 515 #mpa                  
kic = 59.5 #mpa*sqrt(m)          #valor obtido via Granta Design LTDA
c =  6.89*10**(-9) #(MPa*sqrt(m))**m  #assumimos ser perlítico com base em trabalhos de congressos
m1 = 3     #assumimos ser perlítico com base em trabalhos de congressos
 
#momento
def m(x):
    return (w*x/2)*(x - l)
 
#graficos 
class graph():
      def momento():
            plt.style.use("cyberpunk")
            x = arange(0, l + 0.1, 0.1)
            plt.plot(x, m(x))
            plt.xlabel('Comprimento [m]')
            plt.ylabel('Momento Fletor [N.m]')
            mplcyberpunk.add_glow_effects()
            return plt.show()
 
#tensões 
 
r = D/2
 
def σxx_p(p):
      return (p*10**3)*r/(2*t)  #Pa
 
def σxx_m(x):
      d = D - 2*t
      I = (pi/64)*( (D**4) - (d**4) )
      y = r
      return -m(x)*y/I   #Pa
 
#fatores
k = [None]*6 #a,b,c,d,e,f
 
k.insert(0, 4.51*sut**(-0.265)) #fator de superficie #lamidado a frio
k.insert(1, 0.5) #fator de tamanho   #para grandes dimensões kb = 0.5
k.insert(2, 1) #fator de carregamento #para carregamento axial + flexão
k.insert(3, 1) #fator de temperatura 
k.insert(4, 0.753) #fator de confiabilidade #considerando 99.9
k.insert(5, 1) #fator de efeitos diversos
 
kf_flexao = 1 
kf_axial = 1 
 

#limite de endurança
def se():
      se1 = 0.5*sut #limite de endurança
      return k[0]*k[1]*k[2]*k[3]*k[4]*se1
 
#σa' (equivalente de Von Mises)
def σa():
 
      σa_flexao = 0
      σa_axial = (1/2) * ( σxx_p(p = p[1]) - σxx_p(p = p[0] ) )
      
      return ( sqrt( (kf_flexao*σa_flexao + ( kf_axial*σa_axial )/0.85 )**2 ) )/10**6   #MPa
 
#σm' (equivalente de Von Mises)
def σm():
 
      σm_flexao = σxx_m(x = l/2)
      σm_axial = (1/2)* ( σxx_p(p = p[1]) + σxx_p(p = p[0]) )
 
      return sqrt( (kf_flexao*σm_flexao +  kf_axial*σm_axial  )**2 )/10**6         #MPa

  
def σmax():
  return σm() + σa()
 
#criterios de falha
 
class falha_fadiga:
      def solderberg():
            n  = symbols('n')
            n = eq( (( (σa()/se()) + (σm()/sy) )*n) - 1, n )
            return n[0]
      
      def goodman_modificado():
            n  = symbols('n')
            n = eq( (( (σa()/se()) + (σm()/sut) )*n) - 1, n)
            return n[0]
            
      def gerber():
            n  = symbols('n')
            n = eq( (( n*σa()/se() ) + (n*σm()/sut)**2 ) - 1, n )
            return n[1]
      
 
#ciclos
def f():
      return 1.12 

 
sg = σmax()
ac = (1/pi) * ((kic/ (f()*sg) )**2)
af = ac #mm

Δσ = ( σxx_p(p=p[1]) - σxx_p(p=p[0]) )/10**6

cte = 1 - m1/2

N = ( (af**cte) - (ai**cte) ) / ( c*cte*( f()*Δσ*sqrt(pi) )**m1 )
 
 
  #saidas 
print('tensão média {:.4f}MPa, tensão alternada {:.4f}MPa'  .format(σm(), σa()) )
print('Tensão máxima {:.4f}MPa' .format(σmax()))

print('n_solderberg {:.4f}' .format(falha_fadiga.solderberg()) )
print( 'n_goodman {:.4f}' .format(falha_fadiga.goodman_modificado()) )
print( 'n_gerber {:.4f}' .format(falha_fadiga.gerber()))
print('Ciclos {:.4f}M' .format(N/10**6))
graph.momento()
