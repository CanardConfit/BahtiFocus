using CanardConfit.NINA.BahtiFocus.Instructions;
using NINA.Core.Model;
using NINA.Core.Utility;
using NINA.Core.Utility.Notification;
using NINA.Equipment.Equipment.MyCamera;
using NINA.Equipment.Interfaces.Mediator;
using NINA.Equipment.Interfaces.ViewModel;
using NINA.Profile.Interfaces;
using NINA.WPF.Base.Interfaces.Mediator;
using NINA.WPF.Base.ViewModel;
using System;
using System.ComponentModel.Composition;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Input;

namespace CanardConfit.NINA.BahtiFocus.Dockables {
    [Export(typeof(IDockableVM))]
    public class BahtiAnalyserVM : DockableVM, ICameraConsumer {
        
        private readonly IApplicationStatusMediator applicationStatusMediator;
        private readonly ICameraMediator cameraMediator;
        private CancellationTokenSource executeCTS;
        private bool optionsExpanded;
        private ApplicationStatus status;

        [ImportingConstructor]
        public BahtiAnalyserVM(IProfileService profileService, IApplicationStatusMediator applicationStatusMediator, ICameraMediator cameraMediator, IImagingMediator imagingMediator, IFilterWheelMediator fwMediator) : base(profileService) {
            this.applicationStatusMediator = applicationStatusMediator;
            this.cameraMediator = cameraMediator;
            
            Title = "Bahtinov Analyser";
            OptionsExpanded = true;
            
            var dict = new ResourceDictionary();
            dict.Source = new Uri("CanardConfit.NINA.BahtiFocus;component/BahtiVisualImage.xaml", UriKind.RelativeOrAbsolute);
            ImageGeometry = (System.Windows.Media.GeometryGroup) dict["BahtinovSvg"];
            ImageGeometry.Freeze();

            BahtiAnalyser = new BahtiAnalyser(profileService, cameraMediator, imagingMediator, fwMediator);
            
            ExecuteCommand = new AsyncCommand<bool>(
                async () => { 
                    using (executeCTS = new CancellationTokenSource()) {
                        return await Execute(new Progress<ApplicationStatus>(p => Status = p), executeCTS.Token); 
                    }
                },
                _ => BahtiAnalyser.Validate() && cameraMediator.IsFreeToCapture(this));
            
            CancelExecuteCommand = new CommunityToolkit.Mvvm.Input.RelayCommand(() => { try { executeCTS?.Cancel(); } catch (Exception) { } });
            PauseCommand = new CommunityToolkit.Mvvm.Input.RelayCommand(Pause, () => !BahtiAnalyser.IsPausing);
            ResumeCommand = new CommunityToolkit.Mvvm.Input.RelayCommand(Resume);
        }
        
        private ApplicationStatus GetStatus(string stat) {
            return new ApplicationStatus { Source = "BHTA", Status = stat };
        }

        private void Pause() {
            BahtiAnalyser.Pause();
        }

        private void Resume() {
            BahtiAnalyser.Resume();
        }
        
        public async Task<bool> Execute(IProgress<ApplicationStatus> externalProgress, CancellationToken token) {
            try {
                OptionsExpanded = false;
                cameraMediator.RegisterCaptureBlock(this);
                using (var localCTS = CancellationTokenSource.CreateLinkedTokenSource(token)) {
                    await BahtiAnalyser.Execute(externalProgress, localCTS.Token);
                }
            } catch (OperationCanceledException) {
            } catch (Exception ex) {
                Logger.Error(ex);
                Notification.ShowError(ex.Message);
            } finally {
                OptionsExpanded = true;
                cameraMediator.ReleaseCaptureBlock(this);
                externalProgress?.Report(GetStatus(string.Empty));
            }
            return false;
        }
        
        public ApplicationStatus Status {
            get {
                return status;
            }
            set {
                status = value;
                if (string.IsNullOrWhiteSpace(status.Source)) {
                    status.Source = "BHTA";
                }

                RaisePropertyChanged();

                applicationStatusMediator.StatusUpdate(status);
            }
        }

        public IAsyncCommand ExecuteCommand { get; }
        
        public ICommand CancelExecuteCommand { get; }
        
        public ICommand PauseCommand { get; }
        
        public ICommand ResumeCommand { get; }
        
        public BahtiAnalyser BahtiAnalyser { get; }

        public bool OptionsExpanded {
            get => optionsExpanded;
            set {
                optionsExpanded = value;
                RaisePropertyChanged();
            }
        }

        public void UpdateDeviceInfo(CameraInfo deviceInfo) {
            
        }

        public void Dispose() {
        }
    }
}
