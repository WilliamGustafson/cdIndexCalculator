This is a program to compute the [cd-index](https://arxiv.org/pdf/1901.04939.pdf) of a given ranked poset. The input poset
is described by specifying the cover relations. There are also a few provided
modules to calculate particular posets.

Usage: cdIndexCalculator.py [-m MODULE] INPUT

When no module name is provided INPUT should either be a file containing
or be itself a poset specified by listing the
cover relations as follows. Proceeding up ranks, starting at rank 1 and
ending with corank 2, for each element e list the labels of those elements e is
covered by seperated by commas and delimit these element lists by semicolons.
Each rank's list is terminated with a semicolon.
For example, the rank 3 Boolean algebra is described as `1,2;1,3;2,3;;`.

When MODULE is either B_n or Q_n INPUT should be a natural number and it specifies
the rank of the Boolean algebra or cube's face lattice to compute.

When MODULE
is uncrossing a lower interval in the [uncrossing poset](https://arxiv.org/pdf/1406.5671.pdf) is computed.
Elements of the uncrossing poset are pairings on a set {1,...,2n} and for this
module INPUT declares a pairing in the form a1,b1,...,an,bn where ai and bi are
paired. For example, providing 1,4,2,5,3,6 for INPUT calculates the cd-index
of the entire uncrossing poset of order 3.

When MODULE is minor_poset the cd-index of the [minor poset](https://arxiv.org/pdf/2205.01200.pdf) of a given
generator-enriched lattice is computed. INPUT should either be a file
containing or itself be a lattice specified by the cover relations in the following
manner. Label the elements of the lattice 0,1,...,n. List for each element
except the minimum those elements that it covers separated by commas. Seperate the
cover lists with semicolons. For example, the rank 3 Boolean algebra is described
as `0;0;0;1,2;1,3;2,3;4,5,6;`. To specify extra generators that are not
join irreducible provide an argument -g followed by a comma separated
list of the generators to add, identified by their labels in the input.


