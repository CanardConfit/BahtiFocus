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
using System.Windows;
using System.Windows.Media.Imaging;

namespace CanardConfit.NINA.BahtiFocus.Instructions {
    public class BahtiAnalyser: INotifyPropertyChanged {
        private readonly IProfileService _profileService;
        private readonly ICameraMediator _cameraMediator;
        private readonly IImagingMediator _imagingMediator;
        private readonly IFilterWheelMediator _fwMediator;
        private PauseTokenSource _pauseTs;
        private bool _isPaused;
        private IRenderedImage _image;
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
        private BahtinovGrabber.BahtinovGrabber.BahtinovCalc _bahtinovCalc;
        private Rectangle _areaZone;

        public struct Rectangle {
            public int Left { get; }
            public int Top { get; }
            public int Width { get; }
            public int Height { get; }

            public Rectangle(int left, int top, int width, int height) {
                Left = left;
                Top = top;
                Width = width;
                Height = height;
            }
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
            
            CameraInfo = this._cameraMediator.GetInfo();
        }

        public bool Validate() {
            
            var i = new List<string>();

            //Camera
            CameraInfo = _cameraMediator.GetInfo();
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
                    progress.Report(new ApplicationStatus() { Status = $"Capturing new image to solve..." });
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

            // Créer une transformation de mise à l'échelle
            var scaleTransform = new System.Windows.Media.ScaleTransform(
                (double)newWidth / source.PixelWidth, 
                (double)newHeight / source.PixelHeight);

            // Appliquer la transformation à la source de l'image
            var transformedBitmap = new TransformedBitmap(source, scaleTransform);

            return transformedBitmap;
        }
        
        public static BitmapSource GetCenterCrop(BitmapSource source, int width, int height)
        {
            if (source == null)
            {
                throw new ArgumentNullException(nameof(source));
            }

            // Calculer la position de départ pour le recadrage
            int x = (source.PixelWidth - width) / 2;
            int y = (source.PixelHeight - height) / 2;

            // Vérifier que les dimensions demandées ne dépassent pas l'image source
            if (x < 0 || y < 0 || width > source.PixelWidth || height > source.PixelHeight)
            {
                throw new ArgumentException("Les dimensions demandées sont hors limites de l'image source.");
            }

            // Créer une zone de recadrage
            var cropRect = new Int32Rect(x, y, width, height);

            // Créer et retourner le bitmap recadré
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
                        Image = image;

                        double diameter = FocalLength / FocalRatio;

                        BitmapSource calcImage = ResizeBitmap(image.Image, Convert.ToInt32(image.Image.Width / 2), Convert.ToInt32(image.Image.Height / 2));

                        calcImage = GetCenterCrop(calcImage, 800, 800);

                        int realWith = calcImage.PixelWidth * 2;
                        int realHeight = calcImage.PixelHeight * 2;
                        
                        // TODO: Check si le rectangle est correct car me parrait très grand / dézoomé
                        AreaZone = new Rectangle((image.Image.PixelWidth - realWith) / 2, (image.Image.PixelHeight - realHeight) / 2, realWith, realHeight);
                        
                        CalcImage = calcImage;
                        
                        BahtinovGrabber.BahtinovGrabber.BahtinovCalc calc =
                            BahtinovGrabber.BahtinovGrabber.CalculateLines(calcImage, ref bahtinovAngles, diameter, FocalLength, PixelSize);

                        BahtinovCalc = calc;
                        
                        // TODO: Maybe wait a little between captures? Maybe in parameters
                        await CoreUtil.Wait(
                            TimeSpan.FromSeconds(2), token,
                            progress, "Wait");

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
            return new ApplicationStatus { Source = "BHTA", Status = stat };
        }
        
        public ApplicationStatus Status { get => _status; set { _status = value; RaisePropertyChanged(); } }

        public IRenderedImage Image { get => _image; set { _image = value; RaisePropertyChanged(); } }
        
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

        public BahtinovGrabber.BahtinovGrabber.BahtinovCalc BahtinovCalc { get => _bahtinovCalc; set { _bahtinovCalc = value; RaisePropertyChanged(); } }

        public BitmapSource CalcImage { get => _calcImage; set { _calcImage = value; RaisePropertyChanged(); } }

        public Rectangle AreaZone { get => _areaZone; set { _areaZone = value; RaisePropertyChanged(); } }

        public event PropertyChangedEventHandler PropertyChanged;
        protected void RaisePropertyChanged([CallerMemberName] string propertyName = null) {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}