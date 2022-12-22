# CI Helium

Configuration interaction (CI) is a method for calculating the electronic structure of atoms and molecules. It is particularly useful for treating multi-electron systems, where the electron-electron interactions cannot be neglected.

## Installation
Please install following dependencies
```
Memoize
Polynomials
ClassicalOrthogonalPolynomials
CGcoefficient
Optim
Infinity
AssociatedLegendrePolynomials
QuadGK
```

## Usage

The following command will optimize the the eta parameter in order to minimize ground energy
```
julia main.jl --mode optim --n12 4
```
where the n12 parameter means `n12=n1+n2`. States for n12=4 are as follows
```julia
states = [
        [1, 1, 0, 1],
        [1, 2, 0, 1],
        [1, 3, 0, 1],
        [2, 2, 0, 1],
        [2, 2, 1, 1]
    ]
```
where a state has the following form [n1, n2, l, sign].

# References
This is the implementation of CI from book of Jarda Zamastil.
```
@book{zamastil_kvantova_2016,
	title = {Kvantová mechanika a elektrodynamika},
	isbn = {978-80-246-3223-0},
	abstract = {Učebnice se věnuje výkladu základních principů a přibližných metod  kvantové mechaniky s důrazem na využití symetrií pro získání řešení  praktických problémů. Je zde ukázáno, že symetrie nejen umožňují  elegantní řešení problémů, které jsme schopni vyřešit přesně, ale též  výrazně zjednodušují řešení problémů, které jsme nuceni řešit přibližně.Dále  je zde podán, s minimem nutného formalismu, úvod do kvantové  elektrodynamiky a jejímu použití pro fyziku nízkých energií.  Systematicky je zde rozebrán vliv relativistických, magnetických a  kvantově-elektrodynamických efektů na atomová spektra.},
	language = {cs},
	publisher = {Charles University in Prague, Karolinum Press},
	author = {Zamastil, Jaroslav and Benda, Jakub},
	month = sep,
	year = {2016},
	note = {Google-Books-ID: 2cl3DQAAQBAJ},
	keywords = {Science / Mechanics / Dynamics},
}
```