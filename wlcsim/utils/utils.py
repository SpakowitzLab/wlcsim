"""A module for stuff that should be part of Python itself by default but
isn't."""


def make_decorator_factory_args_optional(decorator):
    """Decorator for decorator (factories) that allows them to be used as
    normal decorators. In other words, if you have a decorator that takes
    arguments, then
        @make_decorator_factory_args_optional
        def decorator(arg1=blah, **kwargs)
    will allow you to use the decorator as
        @decorator
        def func():
    instead of only as
        @decorator()
        def func():
    so that you don't have to tell your users to stick with just one.

    WARNING: This make_decorator_factory_args_optional only works if the
    decorator that it's decorating cannot ever be called with a single,
    non-keyword, callable argument in factory mode. i.e. if your decorator can
    be called like
        @decorator(lambda x: x*x)
        def func():
    then DO NOT USE make_decorator_factory_args_optional!"""
    def wrapped_decorator(*args, **kwargs):
        if len(args) == 1 and callable(args[0]) and len(kwargs) == 0:
            return decorator()(args[0])
        else:
            return decorator(*args, **kwargs)
    return wrapped_decorator

@make_decorator_factory_args_optional
def well_behaved_decorator(has_params=False, single_callable_param=False):
    """Decorator for decorators to allow them to preserve docstrings, __name__
    and __dict__ of wrapped functions, so long as the decorator is fairly
    simple. Helpful for documentation purposes.

    If your decorator takes arguments, you'll want to specify this via
    has_params.

    If a decorator doesn't modify function attributes or docstring, and if it
    doesn't take a single, non-keyword-only, callable argument, then it is
    eligible to use this.

    Whenever a decorator is passed a single, callable, non-keyword-only
    parameter, it is likely that the decorator is being called in "wrapper" or
    normal mode, as
        @decorate
        def func():
    If your decorator takes arguments, then you would normally have to call it
    in "factory" mode, as
        @decorate()
        def func():
    You poor users likely will be miffed that they have to remember which one
    to use, but make_decorator_factory_args_optional can allow decorators
    that take optional arguments to be used in the "normal" mode. However, the
    only way to tell if a decorator is being called in the former mode is to
    check if its argument list consists of a single, non-keyword, callable
    argument. We use the same criterion here to correctly extract the identity
    of the function to be wrapped.

    If your decorator does take a single, callable, non-keyword-only parameter
    (screw you), then you can pass single_callable_param=True to specify that
    you would always like your decorator to be treated as if it's in factory
    mode when we decide what the function being wrapped is. You gain back the
    ability to use well_behaved_decorator, but you lose the ability to use your
    decorator in "normal" mode, if you had somehow managed to make that work
    while having a single, positional, callable argument....
    """
    # a decorator with no arguments is just a function that takes a single
    # argument (the function to be decorated) and returns the decorated version
    # of the function
    def wrap_no_params(decorator):
        """Test wrap_no_params"""
        def wrapped_decorator(f):
            """Test wrap_no_params.wrapped_decorator"""
            # gives the wrapped functions properties of the original functions
            g = decorator(f)
            g.__name__ = f.__name__
            g.__doc__ = f.__doc__
            g.__dict__.update(f.__dict__)
            return g
        # gives the wrapped decorator properties of the original decorator
        wrapped_decorator.__name__ = decorator.__name__
        wrapped_decorator.__doc__ = decorator.__doc__
        wrapped_decorator.__dict__.update(decorator.__dict__)
        return wrapped_decorator

    # a decorator with parameters is just a function that takes arguments then
    # returns a decorator with no arguments
    def wrap_with_params(decorator):
        """Factory that returns the decorated decorator if that decorator takes
        arguments. Due to Python's nonsense, we have to manually check if the
        decorator to be decorated is being called with a single callable as its
        argument, as this is our only indication that it is in fact being
        called in normal decorator mode, e.g.
        @decorator
        def func(): ...
        instead of in factory mode, e.g.
        @decorator()
        def func(): ...

        This is because in the former case, we must return the decorated func,
        and in the latter case, we must return a decorator that can decorate
        func.

        Because of this, if you make a decorator that takes a single, callable,
        non-kw argument, then screw you.
        """
        def wrapped_decorator(*args, **kwargs):
            """Centralizes the check to see if the decorator was called in
            factory mode."""
            # the decoration was made in regular mode, so no params are
            # actually used even if the decorator can take them
            if len(args) == 1 and callable(args[0]) and len(kwargs) == 0 \
                    and not single_callable_param:
                return wrap_no_params(decorator)(args[0])
            # otherwise the decoration was made in factory mode, so we must
            # call the decorator itself to get a decorator with no args back
            decorator_given_params = decorator(*args, **kwargs)
            return wrap_no_params(decorator_given_params)
        return wrapped_decorator

    # we are a decorator with params (i.e. a decorator factory), so we return a
    # function that takes a single argument, the decorator to be decorated :)
    if has_params:
        return wrap_with_params
    else:
        return wrap_no_params
