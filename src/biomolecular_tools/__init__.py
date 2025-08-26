# ....................{ BIG BEAR                        }....................
# Warn about type hint violations in *OTHER* packages outside your control;
# only raise exceptions from violations in your package under your control.
# Again, at the very top of your "{your_package}.__init__" submodule:

# from beartype import BeartypeConf  # <-- this isn't your fault
# from beartype.claw import beartype_all, beartype_this_package  # <-- you didn't sign up for this

# beartype_this_package()  # <-- raise exceptions in your code
# beartype_all(conf=BeartypeConf(violation_type=UserWarning))  # <-- emit warnings from other code
