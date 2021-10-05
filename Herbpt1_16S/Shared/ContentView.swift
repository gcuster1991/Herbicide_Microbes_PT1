//
//  ContentView.swift
//  Shared
//
//  Created by Gordon Custer on 10/5/21.
//

import SwiftUI

struct ContentView: View {
    @Binding var document: Herbpt1_16SDocument

    var body: some View {
        TextEditor(text: $document.text)
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView(document: .constant(Herbpt1_16SDocument()))
    }
}
