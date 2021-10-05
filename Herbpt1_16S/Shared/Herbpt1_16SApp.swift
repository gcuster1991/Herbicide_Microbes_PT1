//
//  Herbpt1_16SApp.swift
//  Shared
//
//  Created by Gordon Custer on 10/5/21.
//

import SwiftUI

@main
struct Herbpt1_16SApp: App {
    var body: some Scene {
        DocumentGroup(newDocument: Herbpt1_16SDocument()) { file in
            ContentView(document: file.$document)
        }
    }
}
