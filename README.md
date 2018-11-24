# Conservation Laws

[![License](https://licensebuttons.net/l/by-nc-nd/3.0/88x31.png)](https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)

This is a set of tools for numerically solving different kinds of conservation or balance laws with OpenCL. Currently implemented are:

- linear transport with constant and variable coefficients,
- the nonlinear magnetic induction equation, 
- ideal MHD.

The code is generalised in a form that allows to easily implement and switch between different systems of conservation or balance laws. All numerical aspects are encapsulated in OpenCL kernels. The OpenCL host side is abstracted with the help of [MatCL](https://github.com/MuMPlaCL/MatCL), an OpenCL interface for MathWorks Matlab. This provides
the user with an intuitive and easy way of handling and processing input and output data without any intricate knowledge of OpenCL and allows for interactive development. Usage of the OpenCL kernels is not limited to Matlab.

Exemplary testcases can be found in the `examples` folder. Additional comments are provided in the source files.

This is work in progress. If you have any questions or want to contribute feel free to contact us.


## Prerequisites & Setup

To run the examples the following must be installed:

 - OpenCL Driver (CPU or GPU)
 - OpenCL C++ Headers (e.g. provided by the OpenCL vendors SDKs)
 - Mathworks Matlab
 - [MatCL](https://github.com/MuMPlaCL/MatCL)

 For ease of use you can add `MatCL` to the search path of Matlab.

 ## License

This project is licensed under the terms of the Creative Commons [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode) license.


 ## Disclaimer

Product and company names may be trademarks or registered trademarks of their respective holders.
Use of them does not imply any affiliation with or endorsement by them or their affiliates.
Everything is provided as is and without warranty. Use at your own risk!
