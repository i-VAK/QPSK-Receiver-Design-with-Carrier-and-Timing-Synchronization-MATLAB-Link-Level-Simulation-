# 📡 QPSK Receiver with Carrier & Clock Synchronization

MATLAB link-level simulation of a QPSK receiver implementing a full synchronization chain including coarse carrier frequency offset (CFO) estimation, timing recovery using an interpolating Early–Late loop, matched filtering, and residual phase tracking. Performance is validated using Monte-Carlo SER vs Eb/N0 under practical impairments.

This project focuses on practical PHY receiver synchronization design under carrier offset, clock offset, and noise.

---

## 🎯 Objectives

- Generate framed QPSK transmission with antipodal preamble  
- Model carrier offset, clock offset, AWGN, and frame misalignment  
- Implement carrier and timing synchronization algorithms  
- Recover symbols under realistic impairments  
- Compare simulated SER with theoretical QPSK performance  

Detailed design notes and derivations are available in the included report PDF.

---

## 🧠 Receiver Processing Chain

```
Frame Generation
→ RRC Pulse Shaping
→ Channel + Impairments
→ Matched Filter
→ Coarse CFO Estimation (Squared Preamble + DFT Comb)
→ Carrier Correction
→ Timing Recovery (Interpolating Early–Late Loop)
→ Residual Phase Tracking
→ Symbol Detection
→ SER Measurement
```

---

## 🔧 Key Implementations

### Frame & Modulation
- Gray-mapped QPSK constellation  
- Fixed frame structure with preamble + data  
- Antipodal preamble chosen so squaring removes modulation for CFO estimation  

### Carrier Frequency Offset Estimation
- Squared-preamble method to extract CFO tone  
- DFT comb search over allowed offset range  
- Comb spacing chosen to keep residual error within fine tracker limits  
- Segment averaging for noise robustness  

### Matched Filtering
- Root Raised Cosine (RRC) TX/RX filter pair  
- ISI control and SNR improvement  

### Timing Recovery
- Interpolating Early–Late timing error detector  
- Fractional delay interpolation  
- Window-averaged timing error  
- Threshold-based loop correction  
- Tracks ppm-level clock offsets  

### Residual Phase Tracking
- Preamble-based phase slope estimation  
- Decision-directed block phase correction  
- Progressive constellation de-rotation  

---

## 📊 Simulation & Evaluation

- Monte-Carlo simulation across multiple frames  
- Randomized per frame:
  - data symbols  
  - noise realization  
  - carrier offset  
  - clock offset  
  - start-sample misalignment  

**Outputs:**
- SER vs Eb/N0 curves  
- Comparison with theoretical QPSK SER  
- Synchronization behavior diagnostics  

Observed performance matches theory after synchronization with no visible error floor.

---

## 🛠 Requirements

- MATLAB  
- Signal Processing / Communications Toolbox (recommended)

---

## ▶️ How to Run

```matlab
run main.m
```

Adjust parameters inside the script:
- Eb/N0 range  
- carrier offset (ppm)  
- clock offset (ppm)  
- number of Monte-Carlo frames  

---

## 📄 Reference

Full implementation details, equations, and design rationale are documented in:

**advcomms_report.pdf**
