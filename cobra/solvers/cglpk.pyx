# distutils: libraries=glpk
# cython: embedsignature=True

from glpk cimport *
from libc.stdlib cimport malloc, free
from cpython cimport bool

from tempfile import NamedTemporaryFile as _NamedTemporaryFile  # for pickling
from os import unlink as _unlink
from warnings import warn as _warn

__glpk_version__ = str(glp_version())
_SUPPORTS_MILP = True
solver_name = "cglpk"

__doc__ = """Bindings to GLPK

The GNU Linear Programming Kit (GLPK) is released under the GPL.
The source can be downloaded from http://www.gnu.org/software/glpk/

"""

cdef dict ERROR_CODES = {
    GLP_EBADB: "GLP_EBADB",
    GLP_ESING: "GLP_ESING",
    GLP_ECOND: "GLP_ECOND",
    GLP_EBOUND: "GLP_EBOUND",
    GLP_EFAIL: "GLP_EFAIL",
    GLP_EOBJLL: "GLP_EOBJLL",
    GLP_EOBJUL: "GLP_EOBJUL",
    GLP_EITLIM: "GLP_EITLIM",
    GLP_ETMLIM: "GLP_ETMLIM",
    GLP_ENOPFS: "GLP_ENOPFS",
    GLP_ENODFS: "GLP_ENODFS",
    GLP_EROOT: "GLP_EROOT",
    GLP_ESTOP: "GLP_ESTOP",
    GLP_EMIPGAP: "GLP_EMIPGAP",
    GLP_ENOFEAS: "GLP_ENOFEAS",
    GLP_ENOCVG: "GLP_ENOCVG",
    GLP_EINSTAB: "GLP_EINSTAB",
    GLP_EDATA: "GLP_EDATA",
    GLP_ERANGE: "GLP_ERANGE"
}

cdef dict ERROR_MESSAGES = {
    GLP_EBADB: "invalid basis",
    GLP_ESING: "singular matrix",
    GLP_ECOND: "ill-conditioned matrix",
    GLP_EBOUND: "invalid bounds",
    GLP_EFAIL: "solver failed",
    GLP_EOBJLL: "objective lower limit reached",
    GLP_EOBJUL: "objective upper limit reached",
    GLP_EITLIM: "iteration limit exceeded",
    GLP_ETMLIM: "time limit exceeded",
    GLP_ENOPFS: "no primal feasible solution",
    GLP_ENODFS: "no dual feasible solution",
    GLP_EROOT: "root LP optimum not provided",
    GLP_ESTOP: "search terminated by application",
    GLP_EMIPGAP: "relative mip gap tolerance reached",
    GLP_ENOFEAS: "no primal/dual feasible solution",
    GLP_ENOCVG: "no convergence",
    GLP_EINSTAB: "numerical instability",
    GLP_EDATA: "invalid data",
    GLP_ERANGE: "result out of range"
}

cdef dict METHODS = {
    "auto": GLP_DUALP,
    "primal": GLP_PRIMAL,
    "dual": GLP_DUAL
}

cdef check_error(int result):
    if result == 0:
        return
    if result not in ERROR_CODES:
        raise RuntimeError("glp_simplex failed with unknown error code 0x%x" %
                           result)
    raise RuntimeError("glp_simplex failed with error code %s: %s" %
                       (ERROR_CODES[result], ERROR_MESSAGES[result]))


# Do not want to print out to terminal. Even when not verbose, output
# will be redirected through the hook.
glp_term_out(GLP_OFF)


cdef int hook(void *info, const char *s):
    """function to redirect sdout to python stdout"""
    print(s)
    return 1


cdef class GLP:
    cdef glp_prob *glp
    cdef glp_smcp parameters
    cdef glp_iocp integer_parameters

    # cython related allocation/dellocation functions
    def __cinit__(self):
        self.glp = glp_create_prob()
        glp_set_obj_dir(self.glp, GLP_MAX)  # default is maximize
        glp_init_smcp(&self.parameters)
        glp_init_iocp(&self.integer_parameters)

    def __dealloc__(self):
        glp_delete_prob(self.glp)

    def __init__(self, cobra_model=None):
        cdef int bound_type, index, m, n, n_values, i
        cdef glp_prob *glp
        cdef int *c_rows
        cdef int *c_cols
        cdef double *c_values
        
        # initialize parameters
        self.parameters.msg_lev = GLP_MSG_OFF
        glp_term_hook(NULL, NULL)

        if cobra_model is None:
            return
        glp = self.glp
        m = len(cobra_model.metabolites)
        n = len(cobra_model.reactions)
        glp_add_rows(glp, m)
        glp_add_cols(glp, n)

        metabolite_id_to_index = {r.id: i for i, r
                                  in enumerate(cobra_model.metabolites, 1)}

        linear_constraint_rows = []
        linear_constraint_cols = []
        linear_constraint_values = []

        # set metabolite/consraint bounds
        for index, metabolite in enumerate(cobra_model.metabolites, 1):
            b = float(metabolite._bound)
            c = metabolite._constraint_sense
            if c == 'E':
                bound_type = GLP_FX  # Set metabolite to steady state levels
            elif c == 'L':
                bound_type = GLP_UP  # x < 2 <==> x has an upper bound of 2
            elif c == 'G':
                bound_type = GLP_LO  # x > 2 <==> x has a lower bound of 2
            else:
                raise ValueError("unsupported bound type: %s" % c)
            glp_set_row_bnds(glp, index, bound_type, b, b)
        
        # set reaction/varaiable bounds
        for index, reaction in enumerate(cobra_model.reactions, 1):
            if reaction.variable_kind == "integer":
                if reaction.lower_bound == 0 and reaction.upper_bound == 1:
                    glp_set_col_kind(self.glp, index, GLP_BV)  # binary
                else:
                    glp_set_col_kind(self.glp, index, GLP_IV)
            if reaction.lower_bound == reaction.upper_bound:
                bound_type = GLP_FX
            else:
                bound_type = GLP_DB
            glp_set_col_bnds(glp, index, bound_type,
                             float(reaction.lower_bound),
                             float(reaction.upper_bound))
            glp_set_obj_coef(glp, index, float(reaction.objective_coefficient))

            for metabolite, coefficient in reaction._metabolites.iteritems():
                metabolite_index = metabolite_id_to_index[metabolite.id]
                linear_constraint_rows.append(metabolite_index)
                linear_constraint_cols.append(index)
                linear_constraint_values.append(coefficient)
        
        # set constraint marix
        # first copy the python lists to c arrays
        n_values = len(linear_constraint_rows)
        c_cols = <int *> malloc((n_values + 1) * sizeof(int))
        c_rows = <int *> malloc((n_values + 1) * sizeof(int))
        c_values = <double *> malloc((n_values + 1) * sizeof(double))
        if c_rows is NULL or c_rows is NULL or c_values is NULL:
            raise MemoryError()
        for i in range(n_values):
            c_rows[i + 1] = linear_constraint_rows[i]
            c_cols[i + 1] = linear_constraint_cols[i]
            c_values[i + 1] = float(linear_constraint_values[i])
        # actually set the values
        glp_load_matrix(glp, n_values, c_rows, c_cols, c_values)
        # free the c arrays
        free(c_rows)
        free(c_cols)
        free(c_values)

    # problem creation and modification
    def create_problem(cls, cobra_model, objective_sense="maximize"):
        problem = cls(cobra_model)
        problem.set_objective_sense(objective_sense)
        return problem
    create_problem = classmethod(create_problem)  # decorator does not work

    cpdef change_variable_bounds(self, int index, double lower_bound,
                                 double upper_bound):
        cdef int bound_type = GLP_FX if lower_bound == upper_bound else GLP_DB
        assert index >= 0
        glp_set_col_bnds(self.glp, index + 1, bound_type, lower_bound, upper_bound)

    def change_coefficient(self, int met_index, int rxn_index, double value):
        cdef int col_length, i
        cdef int *indexes
        cdef double *values
        # glpk uses 1 indexing
        met_index += 1
        rxn_index += 1
        # we first have to get the old column
        col_length = glp_get_mat_col(self.glp, rxn_index, NULL, NULL)
        indexes = <int *> malloc((col_length + 2) * sizeof(int))
        values = <double *> malloc((col_length + 2) * sizeof(double))
        if indexes == NULL or values == NULL:
            raise MemoryError()
        glp_get_mat_col(self.glp, rxn_index, indexes, values)
        # search for duplicate
        for i in range(col_length):
            # if a duplicate exists replace that value and exit
            if indexes[i + 1] == met_index:
                values[i + 1] = value
                glp_set_mat_col(self.glp, rxn_index, col_length, indexes, values)
                return
        # need to add a new entry
        indexes[col_length + 1] = met_index
        values[col_length + 1] = value
        glp_set_mat_col(self.glp, rxn_index, col_length + 1, indexes, values)
        free(indexes)
        free(values)

    def solve_problem(self, **solver_parameters):
        cdef int result
        cdef glp_smcp parameters = self.parameters
        cdef glp_iocp integer_parameters = self.integer_parameters
        cdef glp_prob *glp = self.glp

        if "quadratic_component" in solver_parameters:
            q = solver_parameters.pop("quadratic_component")
            if q is not None:
                raise ValueError("quadratic component must be None for glpk")

        for key, value in solver_parameters.items():
            self.set_parameter(key, value)

        # suspend the gil to allow multithreading
        # multithreading must occur with DIFFERENT glp objects
        # calling solve_problem on the same object from 2 different
        # threads at the same time will probably cause problems
        # because glpk itself is not thread safe

        #with nogil:  # we can use this if glpk ever gets thread-safe malloc
        # try to solve the problem with the existing basis
        if glp_simplex(glp, &parameters) != 0:
            glp_adv_basis(glp, 0)
            check_error(glp_simplex(glp, &parameters))
        if self.is_mip():
            self.integer_parameters.tm_lim = self.parameters.tm_lim
            self.integer_parameters.msg_lev = self.parameters.msg_lev
            #self.integer_parameters.tol_bnd = self.parameters.tol_bnd
            #self.integer_parameters.tol_piv = self.parameters.tol_piv
            #with nogil:
            check_error(glp_intopt(glp, &integer_parameters))
        return self.get_status()

    def solve(cls, cobra_model, **kwargs):
        problem = cls.create_problem(cobra_model)
        problem.solve_problem(**kwargs)
        solution = problem.format_solution(cobra_model)
        #cobra_model.solution = solution
        #return {"the_problem": problem, "the_solution": solution}
        return solution
    solve = classmethod(solve)

    def get_status(self):
        cdef int result = glp_mip_status(self.glp) if self.is_mip() \
            else glp_get_status(self.glp)
        if result == GLP_OPT:
            return "optimal"
        if result == GLP_FEAS:
            return glp_get_status(self.glp)
        if result == GLP_UNDEF:
            return "undefined"
        if result == GLP_UNBND:
            return "unbounded"
        if result == GLP_NOFEAS:
            return "infeasible"
        return "failed"

    cpdef set_objective_sense(self, objective_sense):
        objective_sense = objective_sense.lower()
        if objective_sense == "maximize":
            glp_set_obj_dir(self.glp, GLP_MAX)
        elif objective_sense == "minimize":
            glp_set_obj_dir(self.glp, GLP_MIN)
        else:
            raise ValueError("%s is not a valid objective sense" % objective_sense)

    cpdef set_parameter(self, parameter_name, value):
        """set a solver parameter"""
        if parameter_name == "objective_sense":
            self.set_objective_sense(value)
        elif parameter_name == "time_limit":
            self.parameters.tm_lim = int(1000 * value)
        elif parameter_name == "tolerance_feasibility":
            self.parameters.tol_bnd = float(value)
            self.parameters.tol_dj = float(value)
        elif parameter_name == "tolerance_markowitz":
            self.parameters.tol_piv = float(value)
        elif parameter_name == "tolerance_integer":
            self.integer_parameters.tol_int = float(value)
        elif parameter_name == "mip_gap" or parameter_name == "MIP_gap":
            self.integer_parameters.mip_gap = float(value)
        elif parameter_name == "verbose":
            if not value:  # suppress all output
                self.parameters.msg_lev = GLP_MSG_OFF
                glp_term_hook(NULL, NULL)
                return
            glp_term_hook(hook, NULL)
            if value == "err":
                self.parameters.msg_lev = GLP_MSG_ERR
            elif value is True or value == "all":
                self.parameters.msg_lev = GLP_MSG_ALL
            elif value == "normal":
                self.parameters.msg_lev = GLP_MSG_ON
        elif parameter_name == "iteration_limit":
            self.parameters.it_lim = value
        elif parameter_name == "lp_method":
            self.parameters.meth = METHODS[value]
        elif parameter_name == "threads":
            _warn("multiple threads not supported")
        elif parameter_name == "MIP_gap_abs":
            _warn("setting aboslute mip gap not supported")
        else:
            raise ValueError("unknown parameter " + str(parameter_name))

    cpdef get_objective_value(self):
        if self.is_mip():
            return glp_mip_obj_val(self.glp)
        return glp_get_obj_val(self.glp)

    cpdef change_variable_objective(self, int index, double value):
        assert index >= 0
        glp_set_obj_coef(self.glp, index + 1, value)

    cpdef is_mip(self):
        return glp_get_num_int(self.glp) > 0

    def format_solution(self, cobra_model):
        cdef int i, m, n
        cdef glp_prob *glp = self.glp
        Solution = cobra_model.solution.__class__
        status = self.get_status()
        if status != "optimal":  # todo handle other possible
            return Solution(None, status=status)
        solution = Solution(self.get_objective_value(), status=status)
        m = glp_get_num_rows(glp)
        n = glp_get_num_cols(glp)
        x = [0] * n
        if self.is_mip():
            for i in range(1, n + 1):
                    x[i - 1] = glp_mip_col_val(glp, i)
            solution.x_dict = {rxn.id: x[i] for i, rxn
                               in enumerate(cobra_model.reactions)}
            solution.x = x
        else:
            for i in range(1, n + 1):
                x[i - 1] = glp_get_col_prim(glp, i)
            solution.x_dict = {rxn.id: x[i] for i, rxn
                               in enumerate(cobra_model.reactions)}
            solution.x = x
            y = [0] * m
            for i in range(1, m + 1):
                y[i - 1] = glp_get_row_dual(glp, i)
            solution.y_dict = {met.id: y[i] for i, met
                               in enumerate(cobra_model.metabolites)}
            solution.y = y
        return solution

    # make serializable and copyable
    def __getstate__(self):
        cdef int result
        cdef char *name
        tempfile = _NamedTemporaryFile(mode="r", delete=False)
        name = tempfile.name
        tempfile.close()
        result = glp_write_prob(self.glp, 0, name)
        assert result == 0
        with open(name, "r") as infile:
            state = infile.read()
        _unlink(name)
        return state

    def __reduce__(self):
        return (GLP, (), self.__getstate__())

    def __setstate__(self, state):
        cdef int result
        cdef char *name = NULL
        with _NamedTemporaryFile(mode="w", delete=False) as tempfile:
            name = tempfile.name
            tempfile.write(state)
        result = glp_read_prob(self.glp, 0, name)
        assert result == 0
        _unlink(name)

    def __copy__(self):
        other = GLP()
        glp_copy_prob(other.glp, self.glp, GLP_ON)
        other.parameters = self.parameters
        other.integer_parameters = self.integer_parameters
        return other

# wrappers for all the functions at the module level
create_problem = GLP.create_problem
def set_objective_sense(lp, objective_sense="maximize"):
    return lp.set_objective_sense(lp, objective_sense=objective_sense)
cpdef change_variable_bounds(lp, int index,
                             double lower_bound, double upper_bound):
    return lp.change_variable_bounds(index, lower_bound, upper_bound)
cpdef change_variable_objective(lp, int index, double value):
    return lp.change_variable_objective(index, value)
cpdef change_coefficient(lp, int met_index, int rxn_index, double value):
    return lp.change_coefficient(met_index, rxn_index, value)
cpdef set_parameter(lp, parameter_name, value):
    return lp.set_parameter(parameter_name, value)
def solve_problem(lp, **kwargs):
    return lp.solve_problem(**kwargs)
cpdef get_status(lp):
    return lp.get_status()
cpdef get_objective_value(lp):
    return lp.get_objective_value()
cpdef format_solution(lp, cobra_model):
    return lp.format_solution(cobra_model)
solve = GLP.solve
