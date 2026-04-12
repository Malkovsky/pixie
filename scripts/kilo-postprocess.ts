import fs from "node:fs";

type JsonObject = Record<string, unknown>;

const MAX_TOOL_OUTPUT_LINES = 12;
const MAX_TOOL_LINES_IN_FINAL_OUTPUT = 80;

function stripAnsi(input: string): string {
  return input.replace(/\x1B\[[0-9;]*[A-Za-z]/g, "");
}

function asObject(value: unknown): JsonObject | null {
  if (typeof value === "object" && value !== null) {
    return value as JsonObject;
  }
  return null;
}

function asText(value: unknown): string {
  return typeof value === "string" ? value : "";
}

function asDisplayText(value: unknown): string {
  if (typeof value === "string") {
    return value;
  }

  if (typeof value === "number" || typeof value === "boolean") {
    return String(value);
  }

  if (value === null || value === undefined) {
    return "";
  }

  try {
    return JSON.stringify(value);
  } catch {
    return "";
  }
}

function normalizeLines(text: string): string[] {
  return stripAnsi(text)
    .replace(/\r/g, "")
    .split("\n")
    .map((line) => line.trimEnd())
    .filter((line) => line.length > 0);
}

function isToolEvent(event: JsonObject): boolean {
  const type = asText(event.type);
  if (type === "tool_use" || type === "tool_start" || type === "tool_end" || type === "tool_result") {
    return true;
  }

  return type === "say" && asText(event.say) === "tool";
}

function eventTimestamp(event: JsonObject): string {
  const raw = event.timestamp;
  if (typeof raw !== "number" || !Number.isFinite(raw)) {
    return "????????????";
  }
  return new Date(raw).toISOString().slice(11, 23);
}

function summarizeEvent(event: JsonObject): string {
  const type = asText(event.type);

  if (type === "welcome") {
    return "session started";
  }

  if (type === "error") {
    const message = asText(event.message) || asText(event.content);
    return message ? `error: ${message}` : "error";
  }

  if (type === "text") {
    const part = asObject(event.part);
    const text = part ? asText(part.text).trim() : "";
    if (text) {
      return `assistant: ${text}`;
    }
    return "assistant text";
  }

  if (type === "tool_start") {
    const tool = asText(event.tool) || asText(event.name) || "unknown";
    return `tool start: ${tool}`;
  }

  if (type === "tool_use") {
    const part = asObject(event.part);
    const tool = (part && asText(part.tool)) || asText(event.tool) || asText(event.name) || "unknown";
    const state = part ? asObject(part.state) : null;
    const status = state ? asText(state.status) : "";
    const input = state ? asObject(state.input) : null;
    const description = input ? asText(input.description).trim() : "";
    const command = input ? asText(input.command).trim() : "";
    const label = description || command;

    if (status && label) {
      return `tool ${tool} ${status}: ${label}`;
    }
    if (status) {
      return `tool ${tool} ${status}`;
    }
    if (label) {
      return `tool ${tool}: ${label}`;
    }
    return `tool ${tool}`;
  }

  if (type === "tool_end" || type === "tool_result") {
    const tool = asText(event.tool) || asText(event.name) || "unknown";
    return `tool end: ${tool}`;
  }

  if (type === "say") {
    const say = asText(event.say);
    const content = stripAnsi(asText(event.content)).trim();

    if (say === "text") {
      return content ? `assistant: ${content}` : "assistant text";
    }

    if (say === "reasoning") {
      const partial = event.partial === true;
      if (partial) {
        return "reasoning (partial)";
      }
      return content ? `reasoning: ${content}` : "reasoning";
    }

    if (say === "tool") {
      return content ? `tool: ${content}` : "tool";
    }

    if (say === "api_req_started") {
      const metadata = asObject(event.metadata);
      const protocol = metadata ? asText(metadata.apiProtocol) : "";
      return protocol ? `api request started: ${protocol}` : "api request started";
    }

    return say ? `say/${say}` : "say";
  }

  if (type) {
    return `event: ${type}`;
  }

  return "event";
}

function extractToolOutputLines(event: JsonObject): string[] {
  if (asText(event.type) !== "tool_use") {
    return [];
  }

  const part = asObject(event.part);
  if (!part) {
    return [];
  }

  const state = asObject(part.state);
  if (!state) {
    return [];
  }

  const metadata = asObject(state.metadata);
  const rawOutput = asDisplayText(state.output) || (metadata ? asDisplayText(metadata.output) : "");
  if (!rawOutput.trim()) {
    return [];
  }

  const outputLines = normalizeLines(rawOutput);
  if (outputLines.length === 0) {
    return [];
  }

  const limited = outputLines.slice(0, MAX_TOOL_OUTPUT_LINES);
  const result = ["tool output:"];
  for (const line of limited) {
    result.push(`  ${line}`);
  }

  if (outputLines.length > MAX_TOOL_OUTPUT_LINES) {
    result.push(`  ... (${outputLines.length - MAX_TOOL_OUTPUT_LINES} more lines)`);
  }

  return result;
}

function extractText(event: JsonObject): string {
  const type = asText(event.type);

  if (type === "text") {
    const part = asObject(event.part);
    return part ? asText(part.text).trim() : "";
  }

  if (type === "say" && asText(event.say) === "text") {
    return stripAnsi(asText(event.content)).trim();
  }

  return "";
}

function main(): void {
  const eventsPath = process.argv[2] ?? "kilo-events.log";
  const outputPath = process.argv[3] ?? "kilo-output.log";
  const readablePath = process.argv[4] ?? "kilo-readable.log";

  const raw = fs.readFileSync(eventsPath, "utf8");
  const lines = raw.split(/\r?\n/);

  const textParts: string[] = [];
  const readable: string[] = [];
  const toolTrace: string[] = [];

  for (const line of lines) {
    if (!line.trim()) {
      continue;
    }

    try {
      const parsed = JSON.parse(line) as unknown;
      const event = asObject(parsed);
      if (!event) {
        continue;
      }

      const text = extractText(event);
      if (text) {
        textParts.push(text);
      }

      const timestamp = eventTimestamp(event);
      const summary = summarizeEvent(event);
      const summaryLine = `[${timestamp}] ${summary}`;
      readable.push(summaryLine);

      const toolOutputLines = extractToolOutputLines(event);
      for (const line of toolOutputLines) {
        readable.push(`[${timestamp}] ${line}`);
      }

      if (isToolEvent(event)) {
        toolTrace.push(summaryLine);
        for (const line of toolOutputLines) {
          toolTrace.push(`[${timestamp}] ${line}`);
        }
      }
    } catch {
      const stripped = stripAnsi(line).trim();
      if (stripped) {
        readable.push(stripped);
      }
    }
  }

  let finalOutput = "";
  if (textParts.length > 0) {
    finalOutput = textParts[textParts.length - 1];
  } else {
    const fallback = stripAnsi(raw)
      .split(/\r?\n/)
      .map((item) => item.trim())
      .filter((item) => item.length > 0);
    finalOutput = fallback.length > 0 ? fallback[fallback.length - 1] : "";
  }

  if (toolTrace.length > 0) {
    const limitedToolTrace = toolTrace.slice(0, MAX_TOOL_LINES_IN_FINAL_OUTPUT);
    if (toolTrace.length > MAX_TOOL_LINES_IN_FINAL_OUTPUT) {
      limitedToolTrace.push(
        `... (${toolTrace.length - MAX_TOOL_LINES_IN_FINAL_OUTPUT} more tool lines)`
      );
    }

    const toolSection = [`Tool calls:`, ...limitedToolTrace].join("\n");
    finalOutput = finalOutput ? `${finalOutput}\n\n${toolSection}` : toolSection;
  }

  fs.writeFileSync(outputPath, `${finalOutput}\n`, "utf8");
  fs.writeFileSync(readablePath, `${readable.join("\n")}\n`, "utf8");
}

main();
