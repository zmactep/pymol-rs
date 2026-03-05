use pymol_plugin::prelude::*;
use pymol_plugin::pymol_plugin;

pymol_plugin! {
    name: "hello",
    version: "0.2.0",
    description: "Example plugin: registers a 'hello' command",
    commands: [HelloCommand],
}

struct HelloCommand;

impl Command for HelloCommand {
    fn name(&self) -> &str {
        "hello"
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let name = args.get_str(0).unwrap_or("World");
        ctx.print(&format!("Greetings from plug-in system, {}!", name));
        Ok(())
    }

    fn help(&self) -> &str {
        "hello <name>\n\n    Greets the user from the plugin system.\n    Example: hello Alice"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::None]
    }
}
