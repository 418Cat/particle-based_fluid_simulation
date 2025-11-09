## TODO

1. Implement proper OOP **¹**
2. Testing **¹**
3. Better function documentation **²**
4. Improve particle gravity performance **²**

## Features:

1. Opening/Saving file (Top status bar with file options ?) to help with bugfixing **²**
2. Better camera movement control & GPU shader visuals (Performance time chart?) **²**

## Bugfixes:

1. Bboxes causing unknown axis-aligned force applied on particle **²**
    - Seems influenced by the number of bounding boxes per axis
        - The greater the number of bboxes per axis, the greater the force on that axis
    - Using velocity:
        - Bug appears even at high particle/domain bounciness (Seems influenced mainly by particle bounciness)
        - Forces pushes particles towards axis origin
    - Using acceleration:
        - Causes instabilities at high particle/domain bounciness
        - Forces pushes particles away from axis origin, opposed effect when using velocity
    - Causes particles to be attracted to one side of the domain when lots of collisions happen
    - Also happens with radial gravity
        - Causes particles to form a disk and rotate with increasing velocity
2. SegFault when closing program **¹**
3. SegFault with particles gravity **¹**

## Tweaks

1. Add stuff to readme
2. Fix color display of velocity & acceleration - **Finished**
3. Average info ui data over multiple frames - **Finished**

## Later steps

1. Improve implementation of [SPH Fluids in Computer Graphics](https://cgl.ethz.ch/Downloads/Publications/Papers/2014/Sol14a/Sol14a.pdf) **³**
    - Test & optimize gradient computation
2. Implement computer shaders **³**
3. Implement static forces computation


> **¹, ², ³** : Expected change order


