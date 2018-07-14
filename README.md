# simd

## Motivation
For a long time I have been implementing some simd libraries.
I finally decided to use the GLM library (https://github.com/g-truc/glm), but there are some things that I prefer to implement  differently:
1) It only has headers, which at first is charming, but I think it is better to define the bodies of the functions in cpp files in order to be compiled. GLM argues that it has many files precisely to save compiling time by having it all in headers files, but I think it can be avoided by using "extern" in templates.
2) In spite of being centered in glsl, it is true that it allows to easily change to other systems like hlsl, but also, I would prefer that the names of the classes, also respect other system nomenclature.
3) In my opinion, some classes are missing, and others that exist are like "complements" instead of remain as a core classes.

So, this is the beginning of a new implementation keeping the License of GLM
