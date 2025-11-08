## TODO

0. Add stuff to readme

1. Implement proper OOP **¹**
2. Testing **¹**
3. Better funtion documentation **²**
4. Improve particle gravity performance **²**

## Features:

1. Opening/Saving file (Stop satus bar with file options ?) to help with bugfixing **²**
2. Better camera movement control & GPU shader visuals (Performance time chart?) **²**

## Bugfixes:

1. Bboxes causing unknown axis-aligned force applied on particle **²**
    - Seems influenced by the number of bounding boxes per axis
    - Seems to be applied negatively when using acceleration vs. velocity to compute collision physics
    - Causes particles to be attracted to one side of the domain when lots of collisions happen
    - Also happens with radial gravity
        - Causes particles to form a disk and rotate with increasing velocity (p_bounciness ~0.25)

2. SegFault at closing program **¹**

- Improve implementation of [SPH Fluids in Computer Graphics](https://cgl.ethz.ch/Downloads/Publications/Papers/2014/Sol14a/Sol14a.pdf) **³**
    1. Test & optimize gradient computation
- Implement computer shaders **³**


> **¹, ², ³** : Expected change order


