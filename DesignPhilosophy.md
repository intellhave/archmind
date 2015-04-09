# Archmind's Design Philosophy #

Developers and Researchers often find themselves in need of prototyping mesh algorithms.
Therefore, it is crucial to use tools that enable fast and intuitive definitions and easy modifications.
Nevertheless, without a careful design of such a tool, there is the risk of becoming bloated, inefficient or even of posing artificial limits to the user.

Therefore, in designing the archmind API the following requirements were considered:
  * The API should only provide a minimal sub-set of low level construction mesh operators
  * Advanced construction tools can be created by combining low level operators
  * An easy to use interface in python and c++

# Extensibility #

A major consequence of this philosophy is that only a small subset of algorithms and operators will ever be supported by this project.
Therefore, new tools and algorithms will be based on the extensibility capabilities of the project instead of their integration into the core kernel.