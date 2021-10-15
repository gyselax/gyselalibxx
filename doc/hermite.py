import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

n = 5
d = 2

nodes = np.linspace(-d,n+d,n+1+2*d)
print(nodes)

t_sym = sp.Symbol('t')

def f1_scalar(t,n1):
    if t>nodes[n1] and t<=nodes[n1+1]:
        return 1
    else:
        return 0

def f1_sympy(t,n1):
    if n1<d:
        return 0
    elif n1>=(n+d):
        return 0
    else:
        return sp.Heaviside(t-nodes[n1])*(1-sp.Heaviside(t-nodes[n1+1]))

f1 = np.vectorize(f1_scalar)
f1_sp = np.vectorize(f1_sympy)

def f(t,n1,m):
    if m==0:
        if isinstance(t,sp.Symbol):
            return f1_sp(t,n1)
        else:
            return f1(t,n1)
    else:
        part1 = (t-nodes[n1])/(nodes[n1+m]-nodes[n1])*f(t,n1,m-1) if nodes[n1]!=nodes[n1+m] else 0
        part2 = (nodes[n1+m+1]-t)/(nodes[n1+m+1]-nodes[n1+1]) * f(t,n1+1,m-1) if nodes[n1+m+1] != nodes[n1+1] else 0
        return part1 + part2

x = np.linspace(-d,n+d,(n+2*d)*20)
greville_pts = np.array( [np.sum(nodes[i:i+d])/d for i in range(1,1+n+d)] )
greville_pts = np.around(greville_pts, decimals = 15)

plt.figure(figsize=[12,6])
plt.rc('text', usetex=True)
plt.rc('font', size=20)
plt.plot([-3,-3],[0,1.1],'k',label='knots')
for i in range(-2,n+d+1):
    plt.plot([i,i],[0,1.1],'k')
plt.plot([0,0],[0,1.1],'k', linewidth=4)
plt.plot([n,n],[0,1.1],'k', linewidth=4)
for i in range(n+d):
    plt.plot(x,f(x,i,d),label="$B_"+str(i)+"(x)$",color='C'+str(i),linewidth=4)
plt.plot(greville_pts[0], 0, 'ok', label='greville\n points', markersize=10)
for i in greville_pts[1:]:
    plt.plot(i, 0, 'ok', markersize=10)
plt.xticks(nodes)
plt.title('Hermite Boundary Conditions')
print(plt.xlim())
plt.ylim(-0.05,1.25)

ax = plt.gca()
# Shrink current axis by 10%
box = ax.get_position()
print(box)
ax.set_position([box.x0/2, box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#plt.savefig("/home/emily/Documents/wip-docs-and-articles/spline_bcs/Figs/hermite_extended.pdf")

plt.show()
