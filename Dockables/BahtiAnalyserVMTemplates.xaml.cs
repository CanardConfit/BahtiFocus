using System.ComponentModel.Composition;
using System.Windows;

namespace CanardConfit.NINA.BahtiFocus.Dockables {
    [Export(typeof(ResourceDictionary))]
    public partial class BahtiAnalyserVMTemplates {
        public BahtiAnalyserVMTemplates() {
            InitializeComponent();
        }
    }
}