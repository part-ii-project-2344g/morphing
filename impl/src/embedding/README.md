Performs the spherical embedding step as described by Alexa's paper.

**Note 1.**

My banana didn't want to converge once with a custom sphere. 
There was a flipped triangle that converged to a line segment, but from the wrong side.

My good guess is that was a fixed vertex. Changing fixed vertices every now and then, or setting a smaller initial epsilon
should help deal with that problem.