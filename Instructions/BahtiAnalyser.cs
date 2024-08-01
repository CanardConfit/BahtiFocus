using CanardConfit.NINA.BahtiFocus.Dockables;
using CanardConfit.NINA.BahtiFocus.Locale;
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
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media.Imaging;

namespace CanardConfit.NINA.BahtiFocus.Instructions {
    public class BahtiAnalyser: INotifyPropertyChanged {
        private readonly IProfileService _profileService;
        private readonly ICameraMediator _cameraMediator;
        private readonly IImagingMediator _imagingMediator;
        private readonly IFilterWheelMediator _fwMediator;
        private PauseTokenSource _pauseTs;
        private bool _isPaused;
        private BitmapSource _calcImage;
        private ApplicationStatus _status;
        
        private IList<string> _issues = new List<string>();
        private CameraInfo _cameraInfo;
        private FilterInfo _filter;
        private double _exposureTime;
        private int _gain;
        private int _offset;
        private BinningMode _binning;
        private double _focalLength;
        private double _focalRatio;
        private double _pixelSize;
        private Bahtinov.Bahtinov.BahtinovCalc _bahtinovCalc;
        private Rectangle _areaZone;
        private int _zoom;
        private int _x_offset;
        private int _y_offset;
        private BitmapSource _analysedImage;

        public readonly struct Rectangle(int left, int top, int width, int height) {
            public int Left { get; } = left;
            public int Top { get; } = top;
            public int Width { get; } = width;
            public int Height { get; } = height;
        }
        
        [ImportingConstructor]
        public BahtiAnalyser(IProfileService profileService, ICameraMediator cameraMediator, IImagingMediator imagingMediator, IFilterWheelMediator fwMediator) {
            this._profileService = profileService;
            this._cameraMediator = cameraMediator;
            this._fwMediator = fwMediator;
            this._imagingMediator = imagingMediator;
            
            Status = new ApplicationStatus();
            
            Filter = profileService.ActiveProfile.PlateSolveSettings.Filter;
            Gain = profileService.ActiveProfile.PlateSolveSettings.Gain;
            Offset = -1;
            ExposureTime = profileService.ActiveProfile.PlateSolveSettings.ExposureTime;
            Binning = new BinningMode(profileService.ActiveProfile.PlateSolveSettings.Binning, profileService.ActiveProfile.PlateSolveSettings.Binning);
            PixelSize = profileService.ActiveProfile.CameraSettings.PixelSize;
            FocalLength = profileService.ActiveProfile.TelescopeSettings.FocalLength;
            FocalRatio = profileService.ActiveProfile.TelescopeSettings.FocalRatio;

            Zoom = 500;
            XOffset = 0;
            YOffset = 0;

            ResetZoomCommand = new CommunityToolkit.Mvvm.Input.RelayCommand(() => Zoom = 500);
            ResetXOffsetCommand = new CommunityToolkit.Mvvm.Input.RelayCommand(() => XOffset = 0);
            ResetYOffsetCommand = new CommunityToolkit.Mvvm.Input.RelayCommand(() => YOffset = 0);

            SliderValueChanged = new CommunityToolkit.Mvvm.Input.RelayCommand(UpdateArea);
            
            CameraInfo = _cameraMediator.GetInfo();

            BahtinovCalc = new Bahtinov.Bahtinov.BahtinovCalc();
            BahtinovCalc.Angles1 = 0;
            BahtinovCalc.Angles2 = 0;
            BahtinovCalc.Angles3 = 0;
            BahtinovCalc.MaskAngle = 0;
            BahtinovCalc.AbsoluteFocusError = 0;
            BahtinovCalc.CriticalFocusThreshold = 0;
            BahtinovCalc.FocusError = 0;
        }

        private void UpdateArea() {
            
            AreaZone = new Rectangle(((CalcImage.PixelWidth - Zoom) / 2) + XOffset, ((CalcImage.PixelHeight - Zoom) / 2) + YOffset, Zoom, Zoom);
        }

        public bool Validate() {
            
            var i = new List<string>();

            // Camera
            CameraInfo = _cameraMediator.GetInfo();
            if (!CameraInfo.Connected) {
                //i.Add(Loc.Instance["LblCameraNotConnected"]);
                i.Add(BahtiLoc.Instance["LblError"]);
            } else {
                if (CameraInfo.CanSetGain && Gain > -1 && (Gain < CameraInfo.GainMin || Gain > CameraInfo.GainMax)) {
                    i.Add(string.Format(Loc.Instance["Lbl_SequenceItem_Imaging_TakeExposure_Validation_Gain"], CameraInfo.GainMin, CameraInfo.GainMax, Gain));
                }
                if (CameraInfo.CanSetOffset && Offset > -1 && (Offset < CameraInfo.OffsetMin || Offset > CameraInfo.OffsetMax)) {
                    i.Add(string.Format(Loc.Instance["Lbl_SequenceItem_Imaging_TakeExposure_Validation_Offset"], CameraInfo.OffsetMin, CameraInfo.OffsetMax, Offset));
                }
            }

            // Filter wheel
            if (_filter != null && !_fwMediator.GetInfo().Connected) {
                i.Add(Loc.Instance["LblFilterWheelNotConnected"]);
                i.Add("Either connect the filter wheel or clear the filter selection!");
            }

            Issues = i;
            
            return i.Count == 0;
        }

        public void Pause() {
            if (_pauseTs != null) {
                _pauseTs.IsPaused = true;
                RaisePropertyChanged(nameof(IsPaused));
            }
        }
        public void Resume() {
            if (_pauseTs != null) {
                _pauseTs.IsPaused = false;
                RaisePropertyChanged(nameof(IsPaused));
            }
        }

        private async Task<IRenderedImage> CaptureImage(IProgress<ApplicationStatus> progress, CancellationToken token) {
            IRenderedImage image = null;
            do {
                token.ThrowIfCancellationRequested();
                
                var seq = new CaptureSequence() { Binning = Binning, Gain = Gain, ExposureTime = ExposureTime, Offset = Offset, FilterType = Filter, ImageType = CaptureSequence.ImageTypes.SNAPSHOT };
                
                try {
                    progress.Report(new ApplicationStatus { Status = $"Capturing new image to solve..." });
                    image = await _imagingMediator.CaptureAndPrepareImage(seq, new PrepareImageParameters(true, false), token, progress);
                } catch (Exception ex) { Logger.Error(ex); }

                if (image == null) {
                    await CoreUtil.Wait(TimeSpan.FromSeconds(1), token, progress, "Image capture failed. Retrying...");
                }
            } while (image == null);

            return image;
        }
        
        public static BitmapSource ResizeBitmap(BitmapSource source, int newWidth, int newHeight)
        {
            if (source == null)
            {
                throw new ArgumentNullException(nameof(source));
            }

            var scaleTransform = new System.Windows.Media.ScaleTransform(
                (double)newWidth / source.PixelWidth, 
                (double)newHeight / source.PixelHeight);

            var transformedBitmap = new TransformedBitmap(source, scaleTransform);

            return transformedBitmap;
        }

        public BitmapSource GetCrop(BitmapSource source)
        {
            if (source == null)
            {
                throw new ArgumentNullException(nameof(source));
            }

            var cropRect = new Int32Rect(AreaZone.Left, AreaZone.Top, AreaZone.Width, AreaZone.Height);

            var croppedBitmap = new CroppedBitmap(source, cropRect);
            return croppedBitmap;
        }
        
        public async Task Execute(IProgress<ApplicationStatus> externalProgress, CancellationToken token) {
            try {
                using (var localCTS = CancellationTokenSource.CreateLinkedTokenSource(token)) {
                    _pauseTs = new PauseTokenSource();
                    IProgress<ApplicationStatus> progress = new Progress<ApplicationStatus>(p => {
                        externalProgress?.Report(p);
                    });

                    float[] bahtinovAngles = new float[3];
                    
                    do {
                        await WaitIfPaused(localCTS.Token, progress);

                        IRenderedImage image = await CaptureImage(progress, localCTS.Token);
                        
                        BitmapSource calcImage = ResizeBitmap(image.Image, Convert.ToInt32(image.Image.Width / 2), Convert.ToInt32(image.Image.Height / 2));
                        
                        CalcImage = calcImage;
                        
                        UpdateArea();

                        AnalysedImage = GetCrop(calcImage);
                        
                        double diameter = FocalLength / FocalRatio;
                        
                        Bahtinov.Bahtinov.BahtinovCalc calc =
                            Bahtinov.Bahtinov.CalculateLines(AnalysedImage, ref bahtinovAngles, diameter, FocalLength, PixelSize);

                        BahtinovCalc = calc;
                        
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
                await _pauseTs.Token.WaitWhilePausedAsync(token);
                progress?.Report(GetStatus(string.Empty));
                IsPaused = false;
            }
        }
        
        private ApplicationStatus GetStatus(string stat) {
            return new ApplicationStatus { Source = BahtiAnalyserVM.STATUS_SOURCE, Status = stat };
        }
        
        public ApplicationStatus Status { get => _status; set { _status = value; RaisePropertyChanged(); } }
        
        public bool IsPausing { get => _pauseTs?.IsPaused ?? false; }
        
        public bool IsPaused { get => _isPaused; private set { _isPaused = value; RaisePropertyChanged(); } }

        public IList<string> Issues { get => _issues; set { _issues = value; RaisePropertyChanged(); }}

        public CameraInfo CameraInfo { get => _cameraInfo; private set { _cameraInfo = value; RaisePropertyChanged(); } }
        
        public FilterInfo Filter { get => _filter; set { _filter = value; RaisePropertyChanged(); } }

        public double ExposureTime { get => _exposureTime; set { _exposureTime = value; RaisePropertyChanged(); } }

        public int Gain { get => _gain; set { _gain = value; RaisePropertyChanged(); } }

        public int Offset { get => _offset; set { _offset = value; RaisePropertyChanged(); } }

        public BinningMode Binning { get => _binning; set { _binning = value; RaisePropertyChanged(); } }

        public double FocalLength { get => _focalLength; set { _focalLength = value; RaisePropertyChanged(); } }

        public double FocalRatio { get => _focalRatio; set { _focalRatio = value; RaisePropertyChanged(); } }

        public double PixelSize { get => _pixelSize; set { _pixelSize = value; RaisePropertyChanged(); } }

        public Bahtinov.Bahtinov.BahtinovCalc BahtinovCalc { get => _bahtinovCalc; set { _bahtinovCalc = value; RaisePropertyChanged(); } }

        public BitmapSource CalcImage { get => _calcImage; set { _calcImage = value; RaisePropertyChanged(); RaisePropertyChanged("IsCalcImageAvailable"); } }

        public Rectangle AreaZone { get => _areaZone; set { _areaZone = value; RaisePropertyChanged(); } }

        public int Zoom { get => _zoom; set { _zoom = value; RaisePropertyChanged(); } }

        public int XOffset { get => _x_offset; set { _x_offset = value; RaisePropertyChanged(); } }

        public int YOffset { get => _y_offset; set { _y_offset = value; RaisePropertyChanged(); } }
        public ICommand SliderValueChanged { get; }
        public ICommand ResetZoomCommand { get; }
        public ICommand ResetYOffsetCommand { get; }
        public ICommand ResetXOffsetCommand { get; }

        public BitmapSource AnalysedImage { get => _analysedImage; set { _analysedImage = value; RaisePropertyChanged(); } }

        public bool IsCalcImageAvailable => CalcImage != null;
        
        public event PropertyChangedEventHandler PropertyChanged;
        protected void RaisePropertyChanged([CallerMemberName] string propertyName = null) {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}
