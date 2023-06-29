# HW4

## Discussion

1. Drag for a long time - indicates unstable periodicity, with long intervals between waves.
2. Occasional abnormal waves - indicates the presence of spike artifacts.
3. Motion artifact - based on previous assignments: motion-induced noise in the electrocardiogram (ECG) signal.
4. Is it necessary to perform manual detection?

## ECG Preprocessing

Preprocessing the signal using filters, wavelet transforms, and other methods can help remove noise and interference.

* Discrete wavelet transform (DWT)
* Filter bank and regularity

## Base-Line Shiftting

Baseline noise has a lower frequency compared to the ECG signal itself, and the ECG signal contains rich low-frequency components. Therefore, it is not possible to use a low-pass filter to remove baseline drift.

* Median filtering method
* Wavelet transform method
* Algorithmic moving average filtering method
* Empirical Mode Decomposition (EMD) method
