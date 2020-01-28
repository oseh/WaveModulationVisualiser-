using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Accord.Math;
using NAudio.Wave;
using NAudio.CoreAudioApi;

namespace ScottPlotMicrophoneFFT
{
    public partial class Form1 : Form
    {

        // MICROPHONE ANALYSIS SETTINGS
        private int RATE = 44100; // sample rate of the sound card
        private int BUFFERSIZE = (int)Math.Pow(2, 11        ); // must be a multiple of 2

        enum option {FFT, PM, FM, AM};
        option currentState = option.FFT;
        // prepare class objects
        public BufferedWaveProvider bwp;

        public Form1()
        {
            InitializeComponent();
            SetupGraphLabels();
            StartListeningToMicrophone();
            timerReplot.Enabled = true;
        }

        void AudioDataAvailable(object sender, WaveInEventArgs e)
        {
            bwp.AddSamples(e.Buffer, 0, e.BytesRecorded);
        }

        private void Form1_Load(object sender, EventArgs e)
        {
        }

        public void SetupGraphLabels()
        {
            scottPlotUC1.fig.labelTitle = "Microphone PCM Data";
            scottPlotUC1.fig.labelY = "Amplitude (PCM)";
            scottPlotUC1.fig.labelX = "Time (ms)";
            scottPlotUC1.Redraw();

            scottPlotUC2.fig.labelTitle = "Microphone FFT Data";
            scottPlotUC2.fig.labelY = "Power (raw)";
            scottPlotUC2.fig.labelX = "Frequency (KHz)";
            scottPlotUC2.Redraw();

            scottPlotUC3.fig.labelTitle = "Demodulation";
            scottPlotUC3.fig.labelY = "Power (raw)";
            scottPlotUC3.fig.labelX = "Frequency (KHz)";
            scottPlotUC3.Redraw();
        }

        public void StartListeningToMicrophone(int audioDeviceNumber = 0)
        {
            WaveIn wi = new WaveIn();
            wi.DeviceNumber = audioDeviceNumber;
            wi.WaveFormat = new NAudio.Wave.WaveFormat(RATE, 1);
            wi.BufferMilliseconds = (int)((double)BUFFERSIZE / (double)RATE * 1000.0);
            wi.DataAvailable += new EventHandler<WaveInEventArgs>(AudioDataAvailable);
            bwp = new BufferedWaveProvider(wi.WaveFormat);
            bwp.BufferLength = BUFFERSIZE * 2;
            bwp.DiscardOnBufferOverflow = true;
            try
            {
                wi.StartRecording();
            }
            catch
            {
                string msg = "Could not record from audio device!\n\n";
                msg += "Is your microphone plugged in?\n";
                msg += "Is it set as your default recording device?";
                MessageBox.Show(msg, "ERROR");
            }
        }

        private void Timer_Tick(object sender, EventArgs e)
        {
            // turn off the timer, take as long as we need to plot, then turn the timer back on
            timerReplot.Enabled = false;
            PlotLatestData();
            timerReplot.Enabled = true;
        }

        public int numberOfDraws = 0;
        public bool needsAutoScaling = true;
        public void PlotLatestData()
        {
            // check the incoming microphone audio
            int frameSize = BUFFERSIZE;
            var audioBytes = new byte[frameSize];
            bwp.Read(audioBytes, 0, frameSize);

            // return if there's nothing new to plot
            if (audioBytes.Length == 0)
                return;
            if (audioBytes[frameSize - 2] == 0)
                return;

            // incoming data is 16-bit (2 bytes per audio point)
            int BYTES_PER_POINT = 2;

            // create a (32-bit) int array ready to fill with the 16-bit data
            int graphPointCount = audioBytes.Length / BYTES_PER_POINT;

            // create double arrays to hold the data we will graph
            double[] pcm = new double[graphPointCount];
            double[] fft = new double[graphPointCount];
            System.Numerics.Complex[] fht = new System.Numerics.Complex[graphPointCount];
            double[] fftReal = new double[graphPointCount/2];
            double[] fhtReal = new double[graphPointCount / 2];
            double[] fhtIm = new double[graphPointCount / 2];
            double[] IP = new double[graphPointCount];
            double[] IF = new double[graphPointCount];
            double[] sqareWave = new double[graphPointCount];

            for(int i = 0; i < graphPointCount; i++)
            {
                sqareWave[i] = Math.Sin(0.001*i*2*Math.PI) >= 0 ? 20 : 0;
            }

            // populate Xs and Ys with double data
            for (int i = 0; i < graphPointCount; i++)
            {
                // read the int16 from the two bytes
                Int16 val = BitConverter.ToInt16(audioBytes, i * 2);

                // store the value in Ys as a percent (+/- 100% = 200%)
                pcm[i] = (double)(val) / Math.Pow(2,16) * 200.0;
            }


            // calculate the full FFT
            fft = FFT(pcm);

            fht = FHT(pcm);
            //fhtIm = fht.Im();
            //fhtReal= fht.Re();

            IP = instantaneousPhase(fht);
            IF = instantaneousFrequency(fht);
            // determine horizontal axis units for graphs
            double pcmPointSpacingMs = RATE / 1000;
            double fftMaxFreq = RATE / 2;
            double fftPointSpacingHz = fftMaxFreq / graphPointCount;

            // just keep the real half (the other half imaginary)
            Array.Copy(fft, fftReal, fftReal.Length);
            
            // plot the Xs and Ys for both graphs
            scottPlotUC1.Clear();
            scottPlotUC1.PlotSignal(pcm, pcmPointSpacingMs, Color.Blue);
            switch (currentState)
            {
                case option.FFT:
                    scottPlotUC2.Clear();
                    scottPlotUC2.PlotSignal(fftReal, fftPointSpacingHz, Color.Blue);
                    break;
                case option.PM:
                    scottPlotUC2.Clear();
                    scottPlotUC2.PlotSignal(carrierSignal(pcm), pcmPointSpacingMs, Color.Blue);
                    break;
                case option.FM:
                    scottPlotUC2.Clear();
                    scottPlotUC2.PlotSignal(carrierSignal(pcm), pcmPointSpacingMs, Color.Blue);
                    break;
                case option.AM:
                    scottPlotUC2.Clear();
                    scottPlotUC2.PlotSignal(carrierSignal(pcm), pcmPointSpacingMs, Color.Blue);
                    break;
            }

            // optionally adjust the scale to automatically fit the data
            if (needsAutoScaling)
            {
                scottPlotUC1.AxisAuto();
                scottPlotUC2.AxisAuto();
                //needsAutoScaling = false;
            }

            //scottPlotUC1.PlotSignal(Ys, RATE);

            numberOfDraws += 1;
            lblStatus.Text = $"Analyzed and graphed PCM and FFT data {numberOfDraws} times";

            // this reduces flicker and helps keep the program responsive
            Application.DoEvents(); 

        }

        private void autoScaleToolStripMenuItem_Click(object sender, EventArgs e)
        {
            needsAutoScaling = true;
        }

        private void PMMenuItem_Click (object sender, EventArgs e)
        {
            currentState = option.PM;
        }
        private void FFTMenuItem_Click(object sender, EventArgs e)
        {
            currentState = option.FFT;
        }

        private void FMMenuItem_Click(object sender, EventArgs e)
        {
            currentState = option.FM;
        }
        private void AMMenuItem_Click(object sender, EventArgs e)
        {
            currentState = option.AM;
        }


        private void websiteToolStripMenuItem_Click(object sender, EventArgs e)
        {
            System.Diagnostics.Process.Start("https://github.com/swharden/Csharp-Data-Visualization");
        }

        public double[] FFT(double[] data)
        {
            double[] fft = new double[data.Length];
            System.Numerics.Complex[] fftComplex = new System.Numerics.Complex[data.Length];
            for (int i = 0; i < data.Length; i++)
                fftComplex[i] = new System.Numerics.Complex(data[i], 0.0);
            Accord.Math.FourierTransform.FFT(fftComplex, Accord.Math.FourierTransform.Direction.Forward);
            for (int i = 0; i < data.Length; i++)
                fft[i] = fftComplex[i].Magnitude;
            return fft;
        }
      
        public double[] carrierSignal(double [] instant)
        {
            double W = 0;
            double WF = 0;
            double eps = (1/ (float)RATE);
            double[] carrierSignalArray = new double[instant.Length/2];
            double B = 0;
            double frequency_Mod = 0;
            double F = 1000;
            double time = 0;
            bool flagFM = true;
            bool flagPM = true;
            bool flagAM = true;

            for (int i = 0; i < instant.Length/2; i++)
            {
                time = ((double)i * 0.001);
                //time = (i * eps + numberOfDraws*eps);
                frequency_Mod = instant[i];
                W = 2 * Math.PI * time * F;
                WF = 2 * Math.PI * frequency_Mod;
                B = 100 / WF;


                switch (currentState)
                {
                    case option.PM:
                        if(flagPM)
                        {
                            scottPlotUC2.fig.labelTitle = "Pase modulation Data";
                            scottPlotUC2.fig.labelY = "Power (raw)";
                            scottPlotUC2.fig.labelX = "Time (ms)";
                            scottPlotUC2.Redraw();
                            flagPM = false;
                            flagFM = true;
                            flagAM = true;

                        }

                        carrierSignalArray[i] = Math.Cos(W + Math.Sin(frequency_Mod));
                        break;
                    case option.FM:
                        if(flagFM)
                        {
                            scottPlotUC2.fig.labelTitle = "Frequency modulation Data";
                            scottPlotUC2.fig.labelY = "Power (raw)";
                            scottPlotUC2.fig.labelX = "Time (ms)";
                            scottPlotUC2.Redraw();
                            flagPM = true;
                            flagFM = false;
                            flagAM = true;
                        }

                        carrierSignalArray[i] = Math.Cos(W + Math.Sin(WF));
                        break;
                    case option.AM:
                        if (flagAM)
                        {
                            scottPlotUC2.fig.labelTitle = "Amplitude modulation Data";
                            scottPlotUC2.fig.labelY = "Power (raw)";
                            scottPlotUC2.fig.labelX = "Time (ms)";
                            scottPlotUC2.Redraw();
                            flagPM = true;
                            flagFM = true;
                            flagAM = false;
                        }
                        carrierSignalArray[i] = Math.Cos(W) + frequency_Mod * Math.Cos(WF) * Math.Cos(W);

                        break;
                }
            }

            return carrierSignalArray;
        }

        public double[] instantaneousPhase(System.Numerics.Complex[] complex)
        {
            double[] IP = new double[complex.Length];

            for(int i = 0; i < complex.Length; i++)
            {
                IP[i] = complex[i].Phase;
            }
            return IP;
        }


        public double[] instantaneousFrequency(System.Numerics.Complex[] complex)
        {
            double[] IF = new double[complex.Length];
            double eps = 1 / (float)RATE;
            double dIPdt = 0;
            for (int i = 0; i < complex.Length; i++)
            {
                if(i == complex.Length - 1)
                    dIPdt = (complex[i].Phase) / eps;
                else
                    dIPdt = (complex[i + 1].Phase - complex[i].Phase) / eps;
                IF[i] = (1 / Math.PI * 2) * dIPdt;
            }
            return IF;
        }

        public System.Numerics.Complex[] FHT(double[] data)
        {
            double[] fht = new double[data.Length];
            System.Numerics.Complex[] fhtComplex = new System.Numerics.Complex[data.Length];
            for (int i = 0; i < data.Length; i++)
                fhtComplex[i] = new System.Numerics.Complex(data[i], 0.0);
            Accord.Math.HilbertTransform.FHT(fhtComplex, Accord.Math.FourierTransform.Direction.Forward);
    
            return fhtComplex;
        }

        private void scottPlotUC1_Load(object sender, EventArgs e)
        {

        }
    }
}
