from sympy import Symbol, sympify, simplify, Rational
from copy import deepcopy
import random as rd

# An expression is considered a number if it does not contain variables.
# The function checks whether the given expression is a number.
# Mainly needed to perform all possible actions in reducing a matrix with variables to row-echelon form without the risk of division by zero.
def is_num(expr):
    tmp = simplify(expr)
    for symbol in symbols:
        if symbol in tmp.free_symbols:
            return False
    return True

# Similar to the normalize function in Matrix, it just takes an array of values as input instead of a Matrix.
# The function simplifies expressions in the elements of the array.
def normalize_only_data(data):
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = simplify(data[i][j])

# Permutation class
class Permutation:
    # Initialization of a permutation, takes n and data as input.
    # If data is absent, it initializes the identity permutation of length n.
    def __init__(self, n, data=None):
        if not data:
            self.n = n
            self.a = [i for i in range(n)]
            return

        self.n = n
        self.a = deepcopy(list(data))

        assert (n == len(data))  # Check that the provided n and the length of the permutation match
        assert (len(set(data)) == len(data))  # Check that all elements in the given permutation are distinct
        assert (min(data) == 0 and max(data) == len(data) - 1)  # Check that all elements in the permutation are from 0 to n - 1

    # Copy function
    def __copy__(self):
        return Permutation(self.n, deepcopy(self.a))

    # Function to calculate the sign of the permutation
    def sign(self):
        cnt = 0
        for i in range(self.n):
            for j in range(i + 1, self.n):
                cnt += self.a[i] > self.a[j]
        if cnt % 2 == 0:
            return 1
        else:
            return -1

    # Check for the presence of the next permutation (assuming the permutation of the form {n, n - 1, ..., 1} is the last one.
    # Needed to stop in time during the complete enumeration of all permutations.
    def has_next_permutation(self):
        b = deepcopy(self.a)
        b.sort(reverse=True)
        return self.a != b

    # Function to find the next lexicographic permutation
    def next_permutation(self):
        res = deepcopy(self.a)
        pos = self.n - 1

        while pos > 0 and res[pos - 1] > res[pos]:
            pos -= 1
        pos -= 1

        if pos == -1:
            res = [i for i in range(self.n)]
            return Permutation(self.n, res)

        pos_needed_element = pos + 1
        for i in range(pos + 1, self.n):
            if res[i] > res[pos] and res[pos_needed_element] > res[i]:
                pos_needed_element = i
        res[pos], res[pos_needed_element] = res[pos_needed_element], res[pos]
        res[pos + 1:] = list(sorted(res[pos + 1:]))

        return Permutation(self.n, res)

    # Function to check permutations for equality
    def __eq__(self, other):
        if not isinstance(other, Permutation):  # Check if other is a permutation
            return False

        return self.a == other.a

    # Function for multiplying two permutations
    def __mul__(self, other):
        assert (isinstance(other, Permutation))  # Check if other is a permutation
        assert (self.n == other.n)  # Check if the lengths of the permutations match

        res = Permutation(self.n)
        for i in range(self.n):
            res.a[i] = self.a[other.a[i]]
        return res

    # Auxiliary function for exponentiation of a permutation
    def __binpow__(self, power):
        if power == 0:
            return Permutation(self.n)
        if power % 2:
            return self * self.__binpow__(power - 1)
        else:
            res = self.__binpow__(power // 2)
            return res * res

    # Function for exponentiation of a permutation, called with .pow(n), where n is the exponent
    def __pow__(self, power):
        res = self.__binpow__(power)
        if power < 0:
            res = res.inverse()
        return res

    # Function to get the inverse of the permutation
    def inverse(self):
        res = Permutation(self.n)
        for i in range(self.n):
            res.a[self.a[i]] = i
        return res

    # Function to display the permutation
    def __repr__(self):
        return f'Permutation of {self.n} elements: {[i + 1 for i in self.a]}'

# Matrix class
class Matrix:
    # Initialization of a matrix, takes n, m, data as input.
    # If data is absent, it initializes with zeros.
    def __init__(self, n, m, data=None):
        self.n, self.m = n, m
        self.data = [[0] * m for i in range(n)]

        if isinstance(data, Matrix):
            self.n, self.m = data.n, data.m
            self.data = data.normalize().data
            return

        if data:
            assert (len(data) == n)  # Check if the provided n matches the number of rows in data

            for i in range(n):
                assert (len(
                    data[i]) == m)  # Check if the provided m matches the number of columns in data for all rows

                for j in range(m):
                    self.data[i][j] = simplify(data[i][j])

    # Copy function
    def __copy__(self):
        return Matrix(self.n, self.m, deepcopy(self.data))

    # Function to simplify expressions in the elements of the matrix
    def normalize(self):
        res = deepcopy(self)
        for i in range(self.n):
            for j in range(self.m):
                res.data[i][j] = simplify(res.data[i][j])
        return res

    # Function to create a zero matrix n by m
    def zeroes(n, m):
        return Matrix(n, m, [[0] * m for i in range(n)])

    # Function to create an identity matrix of size n
    def ones(n):
        res = Matrix.zeroes(n, n)
        for i in range(n):
            res.data[i][i] = 1
        return res

    # Function to get the trace of the matrix
    def tr(self):
        assert (self.n == self.m)  # Check if the matrix is square

        res = 0
        for i in range(self.n):
            res += self.data[i][i]
        return simplify(res)

    # Auxiliary function for exponentiation of a matrix
    def __binpow__(self, power):
        assert (self.n == self.m)  # Check if the matrix is square

        if power == 0:
            return Matrix.ones(self.n)
        if power % 2:
            return self * self.__binpow__(power - 1)
        else:
            res = self.__binpow__(power // 2)
            return res * res

    # Function for exponentiation of a matrix, called through .pow(n), where n is the exponent
    def pow(self, power):
        assert (self.n == self.m)  # Check if the matrix is square

        res = self.__binpow__(abs(power))
        if power < 0:
            res = res.inverse()
        return res.normalize()

    # Check for matrix equality
    def __eq__(self, other):
        if not isinstance(other, Matrix):  # Compare with C * E if other is not a matrix
            return self.normalize().data == (other * Matrix.ones(self.n)).normalize().data

        return self.normalize().data == other.normalize().data

    # Function to multiply a matrix by a matrix/number
    # When multiplying a matrix by a number, the number must be on the right side of the matrix
    def __mul__(self, other):
        if isinstance(other, Matrix):
            assert (self.m == other.n)  # Check if multiplication is possible if other is a matrix

            res = Matrix.zeroes(self.n, other.m)
            for i in range(self.n):
                for j in range(other.m):
                    for k in range(self.m):
                        res.data[i][j] += self.data[i][k] * other.data[k][j]
            return res.normalize()

        # Multiplication of a matrix by a number
        C = simplify(other)
        res = deepcopy(self)

        for i in range(self.n):
            for j in range(self.m):
                res.data[i][j] *= C
        return res.normalize()

    # Function to determine the sign
    def __neg__(self):
        return (self * (-1)).normalize()

    # Function to add a matrix to a matrix, or a number to a matrix
    def __add__(self, other):
        if isinstance(other, Matrix):
            assert (
                    self.n == other.n and self.m == other.m)  # Check if the necessary dimensions of the given matrices match

            res = deepcopy(self)
            for i in range(self.n):
                for j in range(self.m):
                    res.data[i][j] += other.data[i][j]
            return res.normalize()

        assert (
                self.n == self.m)  # If a number is passed in other, then add the matrix to C * E, where self must be square

        C = simplify(other) * Matrix.ones(self.n)
        res = deepcopy(self)

        for i in range(self.n):
            for j in range(self.n):
                res.data[i][j] += C.data[i][j]
        return res.normalize()

    # Function to subtract a matrix from a matrix, or a number from a matrix
    def __sub__(self, other):
        return (self + (-other)).normalize()

    # Function to calculate the determinant
    def det(self):
        assert (self.n == self.m)  # Check if the matrix is square

        res = 0

        p = Permutation(self.n)
        while True:
            cur = simplify(p.sign())
            for i in range(self.n):
                cur *= self.data[i][p.a[i]]
            res += cur
            if not p.has_next_permutation():
                break
            p = p.next_permutation()

        return simplify(res)

    # Function to display the matrix in the format required for defining a matrix in typst
    def __repr__(self):
        return 'mat(' + ';\n'.join([', '.join(map(str, i)) for i in self.data]) + ')'

    # Similar to the __repr__ function, but takes an array of values instead of a Matrix
    def show(data):
        return 'mat(' + ';\n'.join([', '.join(map(str, i)) for i in data]) + ')'

    # Auxiliary function for step_view, to display the actions performed and output them in typst format
    def add_action(type, i, j, coef=None):
        res = ' =>^('

        if type == 'div':  # Division by coef
            res += '[' + str(i + 1) + '] = [' + str(j + 1) + '] dot '
            coef = Rational(1, coef)
        else:  # Subtracting one row from another with a coefficient coef
            res += '[' + str(i + 1) + '] -= [' + str(j + 1) + '] dot '

        if not is_num(coef) or coef < 0:  # Print coef in parentheses if it is a negative number or an expression
            res += '(' + str(coef) + ')'
        else:
            res += str(coef)

        res += ') '
        return res

    # Function to get the row-echelon form (REF) of the matrix and the entire process of transformations with the matrix when obtaining its REF in the typst format
    # If use_cols = -1, then when reducing the matrix to REF, all its columns should be used, otherwise only use_cols
    # The function correctly calculates the REF for matrices without variables
    # For matrices containing variables, the function tries to improve their appearance as much as possible, while avoiding possible divisions by zero
    def step_view(self, use_cols=-1, better=True):
        if use_cols < 0:
            use_cols = self.m

        res = deepcopy(self.data)
        typst = Matrix.show(res)
        row, col = 0, 0
        while row < self.n and col < use_cols:
            for i in range(row, self.n):
                if res[i][col] and is_num(res[i][col]):
                    res[i], res[row] = res[row], res[i]
                    break

            if not res[row][col]:
                col += 1
                continue

            if better:
                if is_num(res[row][col]):
                    for j in range(col + 1, self.m):
                        res[row][j] /= res[row][col]
                    typst += Matrix.add_action('div', row, row, res[row][col])
                    res[row][col] /= res[row][col]
                    normalize_only_data(res)
                    typst += Matrix.show(res)

            if not is_num(res[row][col]):
                no_possibility_of_dividing_by_zero = True
                for i in range(row + 1, self.n):
                    if i == row:
                        continue
                    if not is_num(res[row][col]) and not is_num(res[i][col] / res[row][col]):
                        no_possibility_of_dividing_by_zero = False

                if not no_possibility_of_dividing_by_zero:
                    break

                for i in range(self.n):
                    if i == row:
                        continue
                    if is_num(res[i][col] / res[row][col]):
                        coef = res[i][col] / res[row][col]
                        if coef == 0:
                            continue
                        for j in range(col, self.m):
                            res[i][j] -= res[row][j] * coef
                        typst += Matrix.add_action('sub', i, row, coef)
                        normalize_only_data(res)
                        typst += Matrix.show(res)

                row += 1
                col += 1
            else:
                for i in range(self.n):
                    if i == row:
                        continue
                    coef = res[i][col] / res[row][col]
                    if coef == 0:
                        continue
                    for j in range(col, self.m):
                        res[i][j] -= res[row][j] * coef
                    typst += Matrix.add_action('sub', i, row, coef)
                    normalize_only_data(res)
                    typst += Matrix.show(res)
                row += 1
                col += 1
            normalize_only_data(res)

        return Matrix(self.n, self.m, res).normalize(), typst

    # Function to append another matrix to the right of the matrix
    def right_add_matrix(self, other):
        assert (isinstance(other, Matrix))  # Check if other is a matrix
        assert (self.n == other.n)  # Check if the number of rows in matrices matches

        res = deepcopy(self)
        res.m += other.m
        for i in range(self.n):
            for j in range(other.m):
                res.data[i].append(other.data[i][j])

        return res

    # Function to find the inverse of the matrix
    def inverse(self):
        assert (self.n == self.m)  # Check if the matrix is square
        assert (self.det() != 0)  # Check for the existence of the inverse matrix

        AE = self.right_add_matrix(Matrix.ones(self.n))
        AE = AE.step_view()

        res = Matrix.ones(self.n)
        for i in range(self.n):
            for j in range(self.m, AE.m):
                res.data[i][j - self.m] = AE.data[i][j]
        return res

# Function to create a random array of size n by m
def make_random_array(n, m):
    res = [[0] * m for i in range(n)]

    for i in range(n):
        for j in range(m):
            res[i][j] = simplify(rd.randint(1, 10))
            if rd.randint(0, 1) % 2 == 0:
                res[i][j] += a * simplify(rd.randint(1, 10))
            if rd.randint(0, 1) % 2 == 0:
                res[i][j] += b * simplify(rd.randint(1, 10))
    return res

# Function to create a matrix of size n by m with random values, if array is not passed, otherwise - use array
def make_matrix(n, m, array=None):
    a = make_random_array(n, m)
    if array is not None:
        a = array

    mat_a = Matrix(len(a), len(a[0]), a)
    return mat_a
