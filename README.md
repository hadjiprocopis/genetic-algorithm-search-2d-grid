Perl script to solve a problem similar to 
[Highest total sum path problem](https://perlmonks.org/?node_id=11113757), at
[Perlmonks](https://perlmonks.org)

Briefly the problem is as follows:

Cells on a 2D grid contain points to be collected.
Points are integers (negative, positive, zero).
Movement on the grid can be orthogonal. Not diagonal.
In moving from a square to another, one can choose
to collect points from all the squares it visits on route,
or skip the points (the "slide" mode) from all intermediate
squares between source and destination.

Plan a route which aims to maximise the total points collected
(the sum of positive AND negative points of each square en route)
and minimise the distance traveled.

The specific quiz rules are not as important as the general
idea of optimal route planning over 2D.

The main drawback of the Genetic Algorithm approach is that
it does not guarantee the best solution. It just searches
as best as it can but the best route may be missed.

The other minor drawback is that as the grid size increases, the
computational efficiency decreases unavoidably but only slightly
exponentially, it looks more like polynomially.
That's not so bad
given the following indicative timings for a full run of
100 iterations over the specified grid size (and a genetic
population of 500) single-threaded :

16x16 = 1s

32x32 = 3s

64x64 = 6s

128x128 = 18s

256x256 = 40s

512x512 = 85s

1024x1024 = 195s

The number of bits in the chromosome for the 1024x1024 grid was 11265
which is again not bad.

The basic usage is
```./route-planning-2d.pl --grid-width 128 --grid-height 128```

but this will work with various defaults
```./route-planning-2d.pl```


Perl is a great tool!

And [Algorithm::Evolutionary](https://metacpan.org/pod/Algorithm::Evolutionary) by 
[J. J. Merelo-Guervós](https://metacpan.org/author/JMERELO) is very poweful.

bw,

[Andreas Hadjiprocopis / bliako](https://perlmonks.org/?node_id=1165397)
