# Matlab language features used

The run manager code is written in and will be executed within Matlab. The code is organized around Matlab classes, which are part of Matlab's very well developed object oriented programming functionality. The run manager code also makes use of packages, which helps you organize your code so that names don't collide in the Matlab path. If you are familiar with these aspects of Matlab programming, skip ahead to [Key Concepts](concepts.md).

## Matlab Classes

 While there are technical differences, classes in Matlab work similarly to classes in Java and other object-oriented languages, and are very well documented by Mathworks:

 * [Role of Classes in Matlab](https://www.mathworks.com/help/matlab/matlab_oop/classes-in-the-matlab-language.html)
 * [Class Syntax Guide](https://www.mathworks.com/help/matlab/class-syntax-guide.html)

 In writing the small amount of code needed to use the run manager with your own data, you will be overriding a small number of methods in classes that will inherit from the LFADS classes, though **it is not necessary to deeply understand classes in Matlab in order to get up and running**.

If you are not familiar with object-oriented programming, the basic concept is that a class is a sort of fusion between a `struct` type (with its associated data fields, called properties), and a set of associated functions (called methods) that are defined to operate on a class's data. The term _class_ refers to the specification which defines the property names and the methods. A specific variable that holds actual data is referred to as an _instance_, and can be created, manipulated, and passed around in Matlab by its variable name. In the lingo of object-oriented programming, an instance is a variable whose type is some class.

```matlab
y = 3.0;
myInstance = MyClass();
```

In this code, `y` is a normal Matlab variable whose type is `double` (double-precision floating point). `myInstance` is an instance whose type is the class `MyClass`.

For illustration of a complete class definition, consider a class `Multiplier` whose job is simply to multiply numbers by a fixed constant. The definition of the class is located in a file `Mulitplier.m` as follows:

```matlab
classdef Multiplier < handle
    properties
        gain = 1; % constant by which inputs are multiplied
    end

    methods
        % this method is called a constructor, and will be called when creating
        % new instances of this class. Here we provide a way to specify the gain
        % when creating the instance
        function obj = Multiplier(theGain)
            if nargin > 1
                obj.gain = theGain;
            end
        end

        % this method does the actual multiplication. The first argument always refers
        % to the instance variable itself, enabling you to refer to properties and other
        % methods in that instance. Otherwise, the code acts like a normal Matlab function
        function out = multiply(obj, in)
            out = in * obj.gain;
        end
    end
end
```

With this definition complete, we can then use the class at the command line as follows:

```matlab
>> myMult = Multiplier(5);
>> myMult.multiply(10)
50
>> myOtherMult = Multiplier(2);
>> myOtherMult.multiply(10)
20
>> myMult.gain = 3; % only affects myMult, not myOtherMult
>> myMult.multiply(10)
30
>> myOtherMult.multiply(10)
20
```

Here, note that `myMult` is a Matlab variable which holds an instance of the class `Multiplier`. We then assign a value to the property `gain` of this instance, and then call the method `multiply`.

## Matlab Packages

The run manager code is also organized within Matlab packages. Packages are a way of organizing code that are used in many other programming languages, such as Java and Python. In Matlab, a package is simply a folder that begins with a `+`. Within Matlab, you will then refer to these classes by prefixing the class names with the package name. So within Matlab, `LFADS.Run` refers to the class located on the file system at `+LFADS/Run.m`.

The main advantage to using packages is that it keeps the namespace organized. This enables you to have multiple things with the same name on the Matlab path while referring to them uniquely with the package name prefix. To use with LFADS run manager with your data, you will probably want to create your own package to organize this code. So, for example if you had a type of experimental data from a reaching task, you might create a folder somewhere on the Matlab path called `+ReachingTask`, and within it copy the starter code provided from `+LorenzExperiment`. Then you can refer to `ReachingTask.Dataset` and `ReachingTask.Run` from within Matlab.
