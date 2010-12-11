import mswine
import time

# data for 2d deformation
xpo = [[100, 200], [300, 300], [200,100]]         # coordinates
xb = [[100,100], [400,100], [400,400], [100,400]]
xnb = [[50,50], [450,50], [450,450], [50,450]]


def draw_closed_curve(the_canvas, the_x, the_y, fill="black", px_per_segment=100, width=1):
    def k(x):   # curve weight function
        if x!=0:
            return 1/float(x)
        else:
            return 1.0e10  # to avoid zero division

    if len(the_x)<3 or len(the_x)!=len(the_y):
        return
    
    x = the_x+[the_x[i] for i in range(3)]
    y = the_y+[the_y[i] for i in range(3)]
    t = [[i+1] for i in range(len(x))]
    s1 = [[i,i+1] for i in range(1,len(x))]
    fxi = mswine.get_linear_functions(t, x, s1)    # basis functions
    fyi = mswine.get_linear_functions(t, y, s1)    # basis functions
    ox='none'
    oy='none'
    for i in xrange(len(s1)):    # for all simplexes
        if i>0 and i<len(s1)-1: # these simplexes are to grant smoothness only
            for j in xrange(px_per_segment+1):
                ti1 = t[s1[i][0]-1][0]
                ti2 = t[s1[i][1]-1][0]
                ti = ti1 + j*float(ti2 - ti1)/px_per_segment
                F_x=mswine.F_s([ti], t, s1, fxi, k)
                F_y=mswine.F_s([ti], t, s1, fyi, k)
                if ox!='none' and oy!='none':
                    the_canvas.create_line(F_x, F_y, ox, oy, fill=fill, width=width)
                ox = F_x
                oy = F_y


def k(x):   # common weight function
    if x!=0:
        return 1/float(x)
    else:
        return 1.0e10  # to avoid zero division



from Tkinter import *   # initializing graphics
root = Tk()
canvas1 = Canvas(root, height = 512, width = 512, background = "white")

for i in range(len(xb)):
    canvas1.create_line(xb[i][0], xb[i][1], xnb[i][0], xnb[i][1], fill="#880000", arrow="last")

draw_closed_curve(canvas1, [x[0] for x in xpo], [y[1] for y in xpo], px_per_segment=30, fill="#cccccc")

xdb_x=[xnb[i][0]-xb[i][0] for i in range(len(xb))]
xdb_y=[xnb[i][1]-xb[i][1] for i in range(len(xb))]
fxi=mswine.get_constant_functions(xb,xdb_x,[])
fyi=mswine.get_constant_functions(xb,xdb_y,[])
xo=[]
yo=[]
for x in xpo:
    xo+=[x[0]+mswine.F_w(x, xb, [], fxi, k)]
    yo+=[x[1]+mswine.F_w(x, xb, [], fyi, k)]
draw_closed_curve(canvas1, xo, yo, px_per_segment=30)


canvas1.pack({"side": "left"})
root.mainloop()


