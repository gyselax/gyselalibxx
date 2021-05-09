# TODO & co.

## Questions
* how should boundary conditions be handled
  - are they defined by the geometry or operator specific (a bit of both I
    guess)

## TODO
* add info about sliced out dimensions to Block
* add easy classes to hold set of values (coords typically) Seq<...>
* make moving from 1-elem array to scalar easy for coords & co.
* make applying member func on array easy for coords & co.
* support itering over domains
* support exposing sub-views from Block
* Let blocks start at non-zero origin
* Handle periodic dimensions
* Replace use of constants 0, 1, ... to identify dimensions by named constants
* Add singleton to each dimension
* Add versions taking singleton as parameter to operators
