restart
R = QQ[x_0..x_3]
C = minors(2,matrix {{x_0, x_1, x_2}, {x_1, x_2, x_3}})
v = ideal(x_0, x_1^2, x_2^2, x_3)
I = intersect(C,v)
isPrime I --false

isPrime eliminate(x_3, I) -- if we project away from x_3, the projection is still not prime.

randomGLaction =  flatten entries(random(ZZ^4, ZZ^4)*(transpose vars R))
gI = substitute(I, for i from 0 to 3 list x_i => randomGLaction_i) 

isPrime eliminate(x_3, gI) -- true, so a generic projection of the variety defined by I is prime.