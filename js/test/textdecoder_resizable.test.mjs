// Regression guard for the resizable-heap TextDecoder crash.
//
// Emscripten's UTF8ArrayToString decodes the heap IN PLACE: `UTF8Decoder.decode(heapOrArray.subarray(idx,endPtr))`.
// With -sALLOW_MEMORY_GROWTH the WASM heap is a *growable* buffer, and browsers that back it with a RESIZABLE
// ArrayBuffer throw `TextDecoder ... must not be resizable` on that decode — which fires for every embind
// std::string return (e.g. Reader.groupAttrs, read at store OPEN), crashing the viewer. build.sh post-processes
// the emitted glue to copy off the growable/shared heap before decoding. (emsdk 5.x removed the -sTEXTDECODER=0
// manual-loop opt-out, so the flag route is unavailable.)
//
// This test asserts the built artifact carries the guard and no longer contains the raw in-place decode — so a
// toolchain bump or a dropped post-process step (which would silently reintroduce a browser-only open crash that
// no headless test exercises) fails CI here instead. Run AFTER build.sh.
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const dist = join(dirname(fileURLToPath(import.meta.url)), "..", "dist");
const RAW = "UTF8Decoder.decode(heapOrArray.subarray(idx,endPtr))";           // the unguarded, in-place heap decode
const GUARD = /buffer\.resizable\|\|typeof SharedArrayBuffer!="undefined"&&[^)]*instanceof SharedArrayBuffer\)\?heapOrArray\.slice/;

let fail = 0;
for (const m of ["lstar_io.mjs", "lstar_kernels.mjs", "lstar_writer.mjs"]) {
  const src = readFileSync(join(dist, m), "utf8");
  const hasGuard = GUARD.test(src);
  const hasRaw = src.includes(RAW);
  const ok = hasGuard && !hasRaw;
  console.log(`  ${ok ? "ok  " : "FAIL"}  ${m}  guard=${hasGuard} rawInPlaceDecode=${hasRaw}`);
  if (!ok) fail++;
}
console.log(fail === 0 ? "\nPASS — resizable/shared heap is copied off before UTF8 decode in every module"
                       : "\nFAIL — the resizable-heap UTF8 guard is missing (build.sh post-process dropped?)");
process.exit(fail === 0 ? 0 : 1);
