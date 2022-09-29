# [CBI-time-changed Lévy processes for multi-currency modeling](https://doi.org/10.1007/s10479-022-04982-z) in Java


## How to run the examples

Currently we offer two different calibration methodologies:

* Standard calibration via the FFT and the numerical solution of the Riccati equations.
* Deep learning based calibration.

To run the examples, proceed as follows

* The standard calibration can be simply started by using the test class TestCalibration.java
* For the deep calibration you have two choices, you can use the pretrained networks we provide or you can perform your own training of the networks. You need to run TestDataGenerator.java (this can take a while dependeng on the setting you choose) which creates training data, after this you can train the network by running TestDeepApproximation.java. Finally, the calibration is performe by TestDeepCalibration.java





## Dependencies

* [Finmath](https://github.com/finmath/finmath-lib)
* [Deeplearning4j](https://github.com/deeplearning4j/deeplearning4j)


## Reference
[1] Fontana, C. Gnoatto, A., Szulda, G. CBI-time-changed Lévy processes for multi-currency modeling. Ann Oper Res (2022). https://doi.org/10.1007/s10479-022-04982-z
