GEO Library
============

Summary 
-------

Struct and functions to deal with 3D coordinates.

Most important is to compute distance and course between two points
and to compute a new point given a course and distance from another point.

It can handle both Flat Earth and Ellipsoids such as WGS-84 (used by GPS).

Stand-alone functions
---------------------

* bin/geo_example - Example function for using the library.

* bin/geo_test - Test functions to compute using Ellipsoid Earth models.

Build and Run Instructions
------------------

~~~~ 
make example
./bin/geo_example
make test
./bin/geo_test
~~~~ 

Main bibliography
-----------------

* Thaddeus Vincenty. Direct and inverse solutions of geodesics on the ellipsoid with application of nested equations. *Survey review*, 23(176):88–93, 1975. doi: [10.1179/sre.1975.23.176.88](https://doi.org/10.1179/sre.1975.23.176.88).

* Tomás Soler and Larry D Hothem. Important parameters used in geodetic transformations. *Journal of Surveying Engineering*, 115(4):414–417, 1989. doi: [10.1061/(ASCE)0733-9453(1989)115:4(414)](https://doi.org/10.1061/(ASCE)0733-9453(1989)115:4(414)).
