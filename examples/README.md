# ENT examples
*Examples of usage of the ENT library.*

## Description

This folder contains some examples of usage of the ENT library, as:

- [example-physics.c](example-physics.c) showcases how to create physics data
  from DIS structure functions, or PDF data, and how to rescale them to given
  cross-section values.

  > Note that pre-computed physics data are provided for [CMS11][CMS11] and
  > [BGR18][BGR18] DIS cross-sections, that can be directly loaded. Thus, this
  > example is only relevant in the case that different models are needed.

- [example-collide.c](example-collide.c) shows how to generate collisions
  with a target material.

- [example-transport.c](example-transport.c) illustrates how to transport a
  neutrino through a simple geometry, considering a uniform medium of infinite
  extension. An external constrain is applied on the neutrino traveled distance,
  as boundary condition.

## Installation

On UNIX the examples can be compiled from ENT's top [Makefile](../Makefile),
e.g. as:
```bash
make examples
```
The compiled examples are located under the `bin/` folder, e.g. as
`bin/example-physics`.

Note that the examples require ENT physics data in order to work, available as
[GitHub releases][RELEASES] assets.


## License

The examples are provided independently of the ENT library under a separate
public domain license allowing them to be copied without any restriction.


[CMS11]: "https://arxiv.org/abs/1106.3723
[BGR18]: "https://arxiv.org/abs/1808.02034
[RELEASES]: https://github.com/niess/ent/releases
