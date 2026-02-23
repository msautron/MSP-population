class Singleton(type):
    """
    This is a metaclass used to produce a Singleton class, i.e. a class with only one active instance.
    It works as follows.  When a class is instantiated it looks if there is already an instance in the dictionary _instances
    If not it creates the instance and then return it, otherwise it just returns it
    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class RestrictiveSingleton(type):
    """
    Similar to Singleton above, but if an instance already exist retunrs an error
    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(RestrictiveSingleton, cls).__call__(*args, **kwargs)
        else:
            raise RuntimeError(f"Only one instance  of the class {cls.__name__}  can be present at time defined on a single process")
        return cls._instances[cls]