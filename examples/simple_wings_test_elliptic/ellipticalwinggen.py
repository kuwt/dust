import numpy as np

def gen_chordlength(pos,rootchord):
    return rootchord*np.sqrt(1-(pos/half_span)**2)

chord_root = 1
half_span = 4  #b/2
num_pts_half_span = 7
y = np.linspace(0,half_span, num_pts_half_span)
print("y = ", y)
z = gen_chordlength(y, chord_root)
print("z = ", z)

print(gen_chordlength(half_span-0.1,chord_root))
