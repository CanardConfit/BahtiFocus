using Google.Protobuf.WellKnownTypes;
using NINA.Core.Locale;
using NINA.Core.Model;
using NINA.Core.Model.Equipment;
using NINA.Core.Utility;
using NINA.Core.Utility.Notification;
using NINA.Equipment.Equipment.MyCamera;
using NINA.Equipment.Interfaces.Mediator;
using NINA.Equipment.Model;
using NINA.Image.Interfaces;
using NINA.Profile.Interfaces;
using Nito.AsyncEx;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.ComponentModel.Composition;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;

namespace CanardConfit.NINA.BahtiFocus.Instructions {
    public class BahtiAnalyser: INotifyPropertyChanged {
        private readonly IProfileService profileService;
        private readonly ICameraMediator cameraMediator;
        private readonly IImagingMediator imagingMediator;
        private readonly IFilterWheelMediator fwMediator;
        private PauseTokenSource pauseTS;
        private bool isPaused;
        private IRenderedImage image;
        private ApplicationStatus status;
        
        private IList<string> issues = new List<string>();
        private CameraInfo cameraInfo;
        private FilterInfo filter;
        private double exposureTime;
        private int gain;
        private int offset;
        private BinningMode binning;
        private double focalLength;
        private double focalRatio;
        private double pixelSize;
        
        [ImportingConstructor]
        public BahtiAnalyser(IProfileService profileService, ICameraMediator cameraMediator, IImagingMediator imagingMediator, IFilterWheelMediator fwMediator) {
            this.profileService = profileService;
            this.cameraMediator = cameraMediator;
            this.fwMediator = fwMediator;
            this.imagingMediator = imagingMediator;
            
            Status = new ApplicationStatus();
            
            Filter = profileService.ActiveProfile.PlateSolveSettings.Filter;
            Gain = profileService.ActiveProfile.PlateSolveSettings.Gain;
            Offset = -1;
            ExposureTime = profileService.ActiveProfile.PlateSolveSettings.ExposureTime;
            Binning = new BinningMode(profileService.ActiveProfile.PlateSolveSettings.Binning, profileService.ActiveProfile.PlateSolveSettings.Binning);
            PixelSize = profileService.ActiveProfile.CameraSettings.PixelSize;
            FocalLength = profileService.ActiveProfile.TelescopeSettings.FocalLength;
            FocalRatio = profileService.ActiveProfile.TelescopeSettings.FocalRatio;
            
            CameraInfo = this.cameraMediator.GetInfo();
        }

        public bool Validate() {
            
            var i = new List<string>();

            //Camera
            CameraInfo = cameraMediator.GetInfo();
            if (!CameraInfo.Connected) {
                i.Add(Loc.Instance["LblCameraNotConnected"]);
            } else {
                if (CameraInfo.CanSetGain && Gain > -1 && (Gain < CameraInfo.GainMin || Gain > CameraInfo.GainMax)) {
                    i.Add(string.Format(Loc.Instance["Lbl_SequenceItem_Imaging_TakeExposure_Validation_Gain"], CameraInfo.GainMin, CameraInfo.GainMax, Gain));
                }
                if (CameraInfo.CanSetOffset && Offset > -1 && (Offset < CameraInfo.OffsetMin || Offset > CameraInfo.OffsetMax)) {
                    i.Add(string.Format(Loc.Instance["Lbl_SequenceItem_Imaging_TakeExposure_Validation_Offset"], CameraInfo.OffsetMin, CameraInfo.OffsetMax, Offset));
                }
            }

            //Filter wheel
            if (filter != null && !fwMediator.GetInfo().Connected) {
                i.Add(Loc.Instance["LblFilterWheelNotConnected"]);
                i.Add("Either connect the filter wheel or clear the filter selection!");
            }

            Issues = i;
            return i.Count == 0;
        }

        public void Pause() {
            if (pauseTS != null) {
                pauseTS.IsPaused = true;
                RaisePropertyChanged(nameof(IsPaused));
            }
        }
        public void Resume() {
            if (pauseTS != null) {
                pauseTS.IsPaused = false;
                RaisePropertyChanged(nameof(IsPaused));
            }
        }

        private async Task<IRenderedImage> CaptureImage(IProgress<ApplicationStatus> progress, CancellationToken token) {
            IRenderedImage image = null;
            do {
                token.ThrowIfCancellationRequested();
                
                var seq = new CaptureSequence() { Binning = Binning, Gain = Gain, ExposureTime = ExposureTime, Offset = Offset, FilterType = Filter, ImageType = CaptureSequence.ImageTypes.SNAPSHOT };
                
                try {
                    progress.Report(new ApplicationStatus() { Status = $"Capturing new image to solve..." });
                    image = await imagingMediator.CaptureAndPrepareImage(seq, new PrepareImageParameters(true, false), token, progress);
                } catch (Exception ex) { Logger.Error(ex); }

                if (image == null) {
                    await CoreUtil.Wait(TimeSpan.FromSeconds(1), token, progress, "Image capture failed. Retrying...");
                }
            } while (image == null);

            return image;
        }
        
        public async Task Execute(IProgress<ApplicationStatus> externalProgress, CancellationToken token) {
            try {
                using (var localCTS = CancellationTokenSource.CreateLinkedTokenSource(token)) {
                    pauseTS = new PauseTokenSource();
                    IProgress<ApplicationStatus> progress = new Progress<ApplicationStatus>(p => {
                        externalProgress?.Report(p);
                    });

                    //TODO: add initializations here

                    do {
                        await WaitIfPaused(localCTS.Token, progress);

                        IRenderedImage image = await CaptureImage(progress, localCTS.Token);
                        Image = image;

                        // TODO: Compute bahtinov mask values to draw on canvas in BahtiVisualImage.xaml
                        
                        // TODO: Maybe wait a little between captures?
                        await CoreUtil.Wait(
                            TimeSpan.FromSeconds(profileService.ActiveProfile.TelescopeSettings.SettleTime), token,
                            progress, "Settling");

                    } while (!localCTS.IsCancellationRequested);
                }
            } catch (OperationCanceledException) {} catch (Exception ex) {
                Logger.Error(ex);
                Notification.ShowError("Bahtinov Analyser failed - " + ex.Message);
                throw;
            } finally {
                IsPaused = false;
                externalProgress?.Report(GetStatus(string.Empty));
            }
        }
        
        private async Task WaitIfPaused(CancellationToken token, IProgress<ApplicationStatus> progress) {
            if (IsPausing) {
                IsPaused = true;
                progress?.Report(GetStatus("Paused"));
                await pauseTS.Token.WaitWhilePausedAsync(token);
                progress?.Report(GetStatus(string.Empty));
                IsPaused = false;
            }
        }
        
        private ApplicationStatus GetStatus(string stat) {
            return new ApplicationStatus { Source = "BHTA", Status = stat };
        }
        
        public ApplicationStatus Status { get => status; set { status = value; RaisePropertyChanged(); } }

        public IRenderedImage Image { get => image; set { image = value; RaisePropertyChanged(); } }
        
        public bool IsPausing { get => pauseTS?.IsPaused ?? false; }
        
        public bool IsPaused { get => isPaused; private set { isPaused = value; RaisePropertyChanged(); } }

        public IList<string> Issues { get => issues; set { issues = value; RaisePropertyChanged(); }}

        public CameraInfo CameraInfo { get => cameraInfo; private set { cameraInfo = value; RaisePropertyChanged(); } }
        
        public FilterInfo Filter { get => filter; set { filter = value; RaisePropertyChanged(); } }

        public double ExposureTime { get => exposureTime; set { exposureTime = value; RaisePropertyChanged(); } }

        public int Gain { get => gain; set { gain = value; RaisePropertyChanged(); } }

        public int Offset { get => offset; set { offset = value; RaisePropertyChanged(); } }

        public BinningMode Binning { get => binning; set { binning = value; RaisePropertyChanged(); } }

        public double FocalLength { get => focalLength; set { focalLength = value; RaisePropertyChanged(); } }

        public double FocalRatio { get => focalRatio; set { focalRatio = value; RaisePropertyChanged(); } }

        public double PixelSize { get => pixelSize; set { pixelSize = value; RaisePropertyChanged(); } }

        public event PropertyChangedEventHandler PropertyChanged;
        protected void RaisePropertyChanged([CallerMemberName] string propertyName = null) {
            this.PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}