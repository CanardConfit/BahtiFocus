using CanardConfit.NINA.BahtiFocus.Properties;
using NINA.Core.Utility;
using NINA.Plugin;
using NINA.Plugin.Interfaces;
using System.ComponentModel;
using System.ComponentModel.Composition;
using System.Runtime.CompilerServices;

namespace CanardConfit.NINA.BahtiFocus {
    [Export(typeof(IPluginManifest))]
    public class BahtiFocusPlugin : PluginBase, INotifyPropertyChanged {
        [ImportingConstructor]
        public BahtiFocusPlugin() {
            if (Settings.Default.UpdateSettings) {
                Settings.Default.Upgrade();
                Settings.Default.UpdateSettings = false;
                CoreUtil.SaveSettings(Settings.Default);
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;
        protected void RaisePropertyChanged([CallerMemberName] string propertyName = null) {
            this.PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}
