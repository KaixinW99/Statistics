'''
from math import sqrt
#tuple
tup = (1,4,9,16,25,6,7,8,9,10)
#sqrt
tupsqrt = tuple(sqrt(x) for x in tup[0:5])
print(tupsqrt)
#square
tupsqua = {x:x**2 for x in tup[5:]}
print(tupsqua)
for i in tupsqua.keys(): print(str(i)+",", end="")
print()
for c in tupsqua.values(): print(str(c)+",",end="")
print()
'''


'''
import random
uc = 'A B C D E F G H I J K L M N \
      O P Q R S T U V W X Y Z'
uclist = uc.split()                       ###uppercase
lclist = uc.lower().split()               ###lowercase
numlist = [str(i) for i in range(0,10)]   ###number
for _ in range(5):
    plist = uclist+lclist+numlist         ###list for passwords
    ucrand = random.sample(uclist,2)
    lcrand = random.sample(lclist,2)
    numrand = random.sample(numlist,2)
    for item in ucrand+lcrand+numrand:
        plist.remove(item)
    ranrand = random.sample(plist,2)
    frand = ucrand+lcrand+numrand+ranrand
    random.shuffle(frand)                ###reorders the elements of a list (attribute)
    print(''.join(frand))
'''

### ax^2+bx+c=0
from math import sqrt
def poly2fsolver(a,b,c):
    delta = b**2-4*a*c
    if a!=0:
        if delta>0:
            print("It has two different real solutions.")
            return [(-b+sqrt(delta))/2/a,(-b-sqrt(delta))/2/a]
        elif delta == 0:
            print("It has two same real solutions.")
            return [-b/2/a]
        else:
            print("It has two different complex solutions.")
            return [(-b+complex(0,sqrt(-delta)))/2/a,(-b-complex(0,sqrt(-delta)))/2/a]
    elif b!=0:
        print("It is linear problem with one solution.")
        return [-c/b]
    elif c!=0:
        print("There is no solution.")
        return None
    else:
        print("The solutions are all numbers.")
        return None

if __name__ == "__main__":
    sol1=poly2fsolver(1,1,5)
    print(sol1)