using System;
using System.Collections.Generic;
using System.ComponentModel.Composition;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;

namespace CanardConfit.NINA.BahtiFocus {
    [Export(typeof(ResourceDictionary))]
    partial class BahtiVisualImage : ResourceDictionary {
        public BahtiVisualImage() {
            InitializeComponent();
        }
    }
}