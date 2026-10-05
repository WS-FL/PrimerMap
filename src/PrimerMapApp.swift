import AppKit
import Foundation

final class AppDelegate: NSObject, NSApplicationDelegate {
    private var window: NSWindow!
    private var inputField: NSTextField!
    private var outputField: NSTextField!
    private var offlineBox: NSButton!
    private var strictBox: NSButton!
    private var runButton: NSButton!
    private var openButton: NSButton!
    private var logView: NSTextView!
    private var task: Process?

    func applicationDidFinishLaunching(_ notification: Notification) {
        window = NSWindow(contentRect: NSRect(x: 0, y: 0, width: 780, height: 565),
                          styleMask: [.titled, .closable, .miniaturizable, .resizable],
                          backing: .buffered, defer: false)
        window.title = "PrimerMap — Sequencing Primer Design"
        window.center()
        window.minSize = NSSize(width: 690, height: 480)
        guard let content = window.contentView else { return }

        func label(_ text: String, _ frame: NSRect, size: CGFloat = 13) -> NSTextField {
            let x = NSTextField(labelWithString: text)
            x.frame = frame
            x.font = NSFont.systemFont(ofSize: size)
            content.addSubview(x)
            return x
        }
        func button(_ text: String, _ frame: NSRect, _ action: Selector) -> NSButton {
            let b = NSButton(title: text, target: self, action: action)
            b.frame = frame
            content.addSubview(b)
            return b
        }
        if let iconPath = Bundle.main.path(forResource: "PrimerMap-logo-en", ofType: "png"),
           let icon = NSImage(contentsOfFile: iconPath) {
            let iconView = NSImageView(frame: NSRect(x: 27, y: 507, width: 46, height: 46))
            iconView.image = icon
            iconView.imageScaling = .scaleProportionallyUpOrDown
            iconView.setAccessibilityLabel("PrimerMap logo")
            content.addSubview(iconView)
        }
        let title = label("PrimerMap", NSRect(x: 88, y: 526, width: 650, height: 26), size: 20)
        title.font = NSFont.boldSystemFont(ofSize: 20)
        label("Design sequencing primers around annotated sgRNA targets",
              NSRect(x: 88, y: 503, width: 650, height: 22), size: 13)
        label("Input", NSRect(x: 25, y: 474, width: 50, height: 25))
        inputField = NSTextField(frame: NSRect(x: 78, y: 471, width: 565, height: 28))
        inputField.placeholderString = "GenBank file or folder (.gb / .gbk)"
        content.addSubview(inputField)
        button("Browse…", NSRect(x: 653, y: 471, width: 100, height: 29), #selector(chooseInput))
        label("Output", NSRect(x: 25, y: 432, width: 50, height: 25))
        outputField = NSTextField(frame: NSRect(x: 78, y: 429, width: 565, height: 28))
        outputField.stringValue = (FileManager.default.homeDirectoryForCurrentUser.appendingPathComponent("Desktop/PrimerMap_Output")).path
        content.addSubview(outputField)
        button("Browse…", NSRect(x: 653, y: 429, width: 100, height: 29), #selector(chooseOutput))
        offlineBox = NSButton(checkboxWithTitle: "Offline mode (requires complete exons and flanks)", target: nil, action: nil)
        offlineBox.frame = NSRect(x: 78, y: 394, width: 520, height: 25)
        content.addSubview(offlineBox)
        strictBox = NSButton(checkboxWithTitle: "Strict amplicon size: 450–550 bp only", target: nil, action: nil)
        strictBox.frame = NSRect(x: 78, y: 366, width: 520, height: 25)
        content.addSubview(strictBox)
        label("Guide-to-primer gap >100 bp on both sides; primer Tm 58–62 °C.",
              NSRect(x: 78, y: 338, width: 680, height: 23), size: 12)
        label("Otherwise tries up to 700 bp. Missing sequence is verified with NCBI RefSeq; no genome-wide screen.",
              NSRect(x: 78, y: 315, width: 680, height: 23), size: 12)
        runButton = button("Design Primers", NSRect(x: 78, y: 265, width: 145, height: 33), #selector(runJob))
        openButton = button("Open Results", NSRect(x: 237, y: 265, width: 145, height: 33), #selector(openOutput))
        openButton.isEnabled = false
        let scroll = NSScrollView(frame: NSRect(x: 25, y: 25, width: 730, height: 225))
        scroll.hasVerticalScroller = true
        scroll.borderType = .bezelBorder
        logView = NSTextView(frame: scroll.bounds)
        logView.isEditable = false
        logView.font = NSFont.monospacedSystemFont(ofSize: 11, weight: .regular)
        scroll.documentView = logView
        content.addSubview(scroll)
        window.makeKeyAndOrderFront(nil)
        NSApp.activate(ignoringOtherApps: true)
    }

    @objc private func chooseInput() {
        let panel = NSOpenPanel()
        panel.canChooseFiles = true
        panel.canChooseDirectories = true
        panel.allowsMultipleSelection = false
        panel.prompt = "Use"
        if panel.runModal() == .OK, let url = panel.url { inputField.stringValue = url.path }
    }

    @objc private func chooseOutput() {
        let panel = NSOpenPanel()
        panel.canChooseFiles = false
        panel.canChooseDirectories = true
        panel.canCreateDirectories = true
        panel.prompt = "Select Folder"
        if panel.runModal() == .OK, let url = panel.url { outputField.stringValue = url.path }
    }

    private func append(_ text: String) {
        logView.string += text
        logView.scrollToEndOfDocument(nil)
    }

    @objc private func runJob() {
        let input = inputField.stringValue.trimmingCharacters(in: .whitespacesAndNewlines)
        let output = outputField.stringValue.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !input.isEmpty, !output.isEmpty else { append("Choose an input and output location.\n"); return }
        guard FileManager.default.fileExists(atPath: input) else { append("Input path does not exist.\n"); return }
        let engine = Bundle.main.resourceURL!.appendingPathComponent("engine/PrimerMapEngine")
        guard FileManager.default.isExecutableFile(atPath: engine.path) else { append("PrimerMap engine is missing: \(engine.path)\n"); return }
        runButton.isEnabled = false
        openButton.isEnabled = false
        logView.string = "Designing primers locally…\n"
        let process = Process()
        process.executableURL = engine
        var args = ["--input", input, "--output", output]
        if offlineBox.state == .on { args.append("--offline") }
        if strictBox.state == .on { args.append("--strict-500") }
        process.arguments = args
        let pipe = Pipe()
        process.standardOutput = pipe
        process.standardError = pipe
        pipe.fileHandleForReading.readabilityHandler = { [weak self] handle in
            let data = handle.availableData
            if !data.isEmpty, let string = String(data: data, encoding: .utf8) {
                DispatchQueue.main.async { self?.append(string) }
            }
        }
        do {
            try process.run()
            task = process
            DispatchQueue.global(qos: .userInitiated).async { [weak self] in
                process.waitUntilExit()
                DispatchQueue.main.async {
                    pipe.fileHandleForReading.readabilityHandler = nil
                    self?.append("\nFinished (exit code \(process.terminationStatus)). Results saved to: \(output)\n")
                    self?.runButton.isEnabled = true
                    self?.openButton.isEnabled = true
                    self?.task = nil
                }
            }
        } catch {
            append("Could not start PrimerMap: \(error)\n")
            runButton.isEnabled = true
        }
    }

    @objc private func openOutput() {
        NSWorkspace.shared.open(URL(fileURLWithPath: outputField.stringValue))
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool { true }
}

let app = NSApplication.shared
let delegate = AppDelegate()
app.delegate = delegate
app.run()
