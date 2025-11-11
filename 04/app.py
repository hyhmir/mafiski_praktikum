import tkinter as tk
import rust as rs
import sys


# print(sys.argv[1:])
# rs.branje(sys.argv[1:])


import wave
import sys
import time
import numpy as np
import sounddevice as sd
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# -------- helper: load wav into numpy array --------
def load_wav(path):

    rs.branje([path, '04/wav.txt'])

    # convert raw bytes to int16 numpy, shape (nframes, nch)
    data = np.loadtxt('04/wav.txt', dtype=np.int16)

    return data

# -------- main: player + tkinter spectrum --------
def play_and_visualize(wav_path):
    data = load_wav(wav_path)
    total_frames = data.shape[0]
    print(f"Loaded: {wav_path}, frames={total_frames}")

    # if stereo, convert to float32 interleaved representation for sounddevice
    audio_out = (data.astype(np.float32) / 32768.0)

    sr = 44100
    nch = 1
    # start playback non-blocking
    # sd.default.samplerate = sr
    # sd.default.channels = nch

    # sounddevice.play accepts numpy arrays; it will stream them
    stream = sd.play(audio_out, samplerate=sr)

    # Create Tkinter window and matplotlib figure
    root = tk.Tk()
    root.title("Simple WAV player + spectrum (updates every 0.1s)")

    fig = Figure(figsize=(6,3))
    ax = fig.add_subplot(111)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude")
    ax.set_xlim(-sr/2, sr/2)
    ax.set_ylim(0, 1)  # will autoscale below if needed

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)

    # variables for update
    window_seconds = 0.01  # show 0.1s worth of samples (4410 samples at 44100 Hz)
    window_size = int(sr * window_seconds)  # -> 4410 for 44.1kHz
    hop_ms = 10  # ms update interval (0.1 s)

    # initial empty line
    freqs = np.fft.fftfreq(window_size, d=1.0/sr)
    line_plot = ax.plot(freqs, np.zeros_like(freqs))[0]

    start_time = time.time()

    def update():
        # compute elapsed frames from start_time
        elapsed = time.time() - start_time
        start_frame = int(elapsed * sr)
        end_frame = start_frame + window_size

        # ensure bounds
        if start_frame >= total_frames:
            # finished: stop playback and exit GUI
            try:
                sd.stop()
            except Exception:
                pass
            root.quit()
            return

        # slice frames for FFT; if stereo, mixdown to mono for display
        chunk = data[start_frame:end_frame]
        if chunk.size < window_size:
            # zero-pad tail
            chunk = np.pad(chunk, (0, window_size - chunk.size))
        

        # convert to float, apply window
        x = (chunk.astype(np.float32) / 32768.0)

        spec = np.abs(rs.dft(x.tolist()))
        # normalize for plotting (avoid dividing by zero)
        peak = spec.max() if spec.size and spec.max() != 0 else 1.0
        spec = spec / peak

        # update plot
        line_plot.set_ydata(spec)
        # optionally autoscale y
        ax.set_ylim(0, 1.05)
        canvas.draw_idle()

        # schedule next update
        root.after(hop_ms, update)

    # start periodic updates
    root.after(0, update)

    # start Tk event loop (blocks until window closed)
    try:
        root.mainloop()
    finally:
        # ensure stop audio
        try:
            sd.stop()
        except Exception:
            pass

# -------- command line entry --------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python player_with_spectrum.py file.wav")
        sys.exit(1)
    play_and_visualize(sys.argv[1])


