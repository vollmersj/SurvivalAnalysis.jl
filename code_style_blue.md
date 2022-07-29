The coding style in this package tries to stick to [Blue](https://github.com/invenia/BlueStyle).

Naming

* [x] UpperCamelCase for modules and type names
* [x] lower_snake_case for method names
* [x] Interpretable names (without abbreviation) & without being verbose
* [x] No method names that are clear from type

Spacing

* [x] 4 spaces per indentation level, no tabs
* [x] whitespace for readability
* [x] no trailing whitespace
* [x] 92 character line length limit
* [x] No padding brackets with spaces

Comments

* [x] necessary comments to explain code
* [x] "Julia" = language; "julia" = executable

Functions

* [x] Include explicit `return`s
* [x] functions mutating at least one of their arguments end in `!`

Custom styles

* [x] Import modules with `using`, with one module per line and at the top of the file
* [x] Never say 'error'!