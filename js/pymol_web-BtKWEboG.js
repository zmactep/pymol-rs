class T {
  static __wrap(e) {
    e = e >>> 0;
    const t = Object.create(T.prototype);
    return t.__wbg_ptr = e, U.register(t, t.__wbg_ptr, t), t;
  }
  __destroy_into_raw() {
    const e = this.__wbg_ptr;
    return this.__wbg_ptr = 0, U.unregister(this), e;
  }
  free() {
    const e = this.__destroy_into_raw();
    o.__wbg_webviewer_free(e, 0);
  }
  /**
   * Count atoms matching a selection expression.
   * @param {string} selection
   * @returns {number}
   */
  count_atoms(e) {
    try {
      const a = o.__wbindgen_add_to_stack_pointer(-16), s = p(e, o.__wbindgen_export, o.__wbindgen_export2), f = g;
      o.webviewer_count_atoms(a, this.__wbg_ptr, s, f);
      var t = w().getInt32(a + 0, !0), _ = w().getInt32(a + 4, !0), c = w().getInt32(a + 8, !0);
      if (c)
        throw d(_);
      return t >>> 0;
    } finally {
      o.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
   * Create a new WebViewer bound to a `<canvas>` element.
   *
   * This is async because WebGPU initialization requires awaiting the
   * adapter and device.
   * @param {string} canvas_id
   * @returns {Promise<WebViewer>}
   */
  static create(e) {
    const t = p(e, o.__wbindgen_export, o.__wbindgen_export2), _ = g, c = o.webviewer_create(t, _);
    return d(c);
  }
  /**
   * Execute a PyMOL command string. Returns JSON with output messages.
   * @param {string} command
   * @returns {any}
   */
  execute(e) {
    const t = p(e, o.__wbindgen_export, o.__wbindgen_export2), _ = g, c = o.webviewer_execute(this.__wbg_ptr, t, _);
    return d(c);
  }
  /**
   * Get projected screen-space labels for overlay rendering.
   *
   * Returns a JSON array of `{ x, y, text, kind }` where coordinates
   * are in physical pixels (divide by `devicePixelRatio` for CSS pixels).
   * @returns {any}
   */
  get_labels() {
    const e = o.webviewer_get_labels(this.__wbg_ptr);
    return d(e);
  }
  /**
   * Get current movie state as JSON.
   * @returns {any}
   */
  get_movie_state() {
    const e = o.webviewer_get_movie_state(this.__wbg_ptr);
    return d(e);
  }
  /**
   * Get info about a specific object.
   * @param {string} name
   * @returns {any}
   */
  get_object_info(e) {
    const t = p(e, o.__wbindgen_export, o.__wbindgen_export2), _ = g, c = o.webviewer_get_object_info(this.__wbg_ptr, t, _);
    return d(c);
  }
  /**
   * Get loaded object names as a JSON array.
   * @returns {any}
   */
  get_object_names() {
    const e = o.webviewer_get_object_names(this.__wbg_ptr);
    return d(e);
  }
  /**
   * Get all named selections as JSON array.
   * @returns {any}
   */
  get_selection_list() {
    const e = o.webviewer_get_selection_list(this.__wbg_ptr);
    return d(e);
  }
  /**
   * Get sequence data for all loaded molecules as JSON.
   * @returns {any}
   */
  get_sequence_data() {
    const e = o.webviewer_get_sequence_data(this.__wbg_ptr);
    return d(e);
  }
  /**
   * Load molecular or map data from bytes.
   *
   * `format` should be one of: "pdb", "xyz", "cif", "mmcif", "bcif",
   * "ccp4", "map", "mrc", "prs"
   *
   * Gzip-compressed data is automatically detected and decompressed.
   * @param {Uint8Array} data
   * @param {string} name
   * @param {string} format
   */
  load_data(e, t, _) {
    try {
      const s = o.__wbindgen_add_to_stack_pointer(-16), f = me(e, o.__wbindgen_export), l = g, H = p(t, o.__wbindgen_export, o.__wbindgen_export2), N = g, X = p(_, o.__wbindgen_export, o.__wbindgen_export2), Y = g;
      o.webviewer_load_data(s, this.__wbg_ptr, f, l, H, N, X, Y);
      var c = w().getInt32(s + 0, !0), a = w().getInt32(s + 4, !0);
      if (a)
        throw d(c);
    } finally {
      o.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
   * Returns true when the scene has changed and needs a re-render.
   * @returns {boolean}
   */
  needs_redraw() {
    return o.webviewer_needs_redraw(this.__wbg_ptr) !== 0;
  }
  /**
   * @param {number} x
   * @param {number} y
   * @param {number} button
   * @param {number} modifiers
   */
  on_mouse_down(e, t, _, c) {
    o.webviewer_on_mouse_down(this.__wbg_ptr, e, t, _, c);
  }
  /**
   * @param {number} x
   * @param {number} y
   * @param {number} modifiers
   */
  on_mouse_move(e, t, _) {
    o.webviewer_on_mouse_move(this.__wbg_ptr, e, t, _);
  }
  /**
   * @param {number} x
   * @param {number} y
   * @param {number} button
   */
  on_mouse_up(e, t, _) {
    o.webviewer_on_mouse_up(this.__wbg_ptr, e, t, _);
  }
  /**
   * @param {number} delta_y
   * @param {number} modifiers
   */
  on_wheel(e, t) {
    o.webviewer_on_wheel(this.__wbg_ptr, e, t);
  }
  /**
   * Ray-cast at physical-pixel canvas coordinates and, when picking is
   * enabled, update the `sele` named selection accordingly.
   *
   * `screen_x` and `screen_y` must be in **physical pixels** (i.e.
   * CSS offset coordinates multiplied by `devicePixelRatio`).
   *
   * Returns a JSON `PickHitInfo` object on a hit, or `null` on a miss or
   * when picking is disabled.
   * @param {number} screen_x
   * @param {number} screen_y
   * @returns {any}
   */
  pick_at_screen(e, t) {
    const _ = o.webviewer_pick_at_screen(this.__wbg_ptr, e, t);
    return d(_);
  }
  /**
   * Update hover indicators by ray-casting at physical-pixel coordinates.
   *
   * Call this from `mousemove` (coordinates pre-multiplied by `devicePixelRatio`).
   * Pass `(-1, -1)` on `mouseleave` — the ray will miss everything and all
   * hover indicators are cleared.
   *
   * No-ops while any mouse button is held (camera drag in progress).
   * @param {number} screen_x
   * @param {number} screen_y
   */
  process_hover(e, t) {
    o.webviewer_process_hover(this.__wbg_ptr, e, t);
  }
  /**
   * Process accumulated input deltas and update the camera.
   *
   * Call this once per frame before `render_frame()`.
   */
  process_input() {
    o.webviewer_process_input(this.__wbg_ptr);
  }
  /**
   * Render one frame to the canvas.
   */
  render_frame() {
    o.webviewer_render_frame(this.__wbg_ptr);
  }
  /**
   * Handle canvas resize.
   * @param {number} width
   * @param {number} height
   */
  resize(e, t) {
    o.webviewer_resize(this.__wbg_ptr, e, t);
  }
  /**
   * Enable or disable cursor-based atom picking (default: disabled).
   *
   * When enabled, `pick_at_screen` performs ray-cast picking and updates
   * the `sele` named selection on click.
   * @param {boolean} enabled
   */
  set_picking_enabled(e) {
    o.webviewer_set_picking_enabled(this.__wbg_ptr, e);
  }
  /**
   * Advance movie playback, rock animation, and camera interpolation.
   *
   * Call this once per frame before `needs_redraw()`.
   * `dt` is the elapsed time in seconds since the last frame.
   * @param {number} dt
   */
  update_animations(e) {
    o.webviewer_update_animations(this.__wbg_ptr, e);
  }
}
Symbol.dispose && (T.prototype[Symbol.dispose] = T.prototype.free);
function he() {
  o.init();
}
function E() {
  return {
    __proto__: null,
    "./pymol_web_bg.js": {
      __proto__: null,
      __wbg_Error_83742b46f01ce22d: function(e, t) {
        const _ = Error(u(e, t));
        return i(_);
      },
      __wbg_Window_e0df001eddf1d3fa: function(e) {
        const t = n(e).Window;
        return i(t);
      },
      __wbg_WorkerGlobalScope_d731e9136c6c49a0: function(e) {
        const t = n(e).WorkerGlobalScope;
        return i(t);
      },
      __wbg___wbindgen_debug_string_5398f5bb970e0daa: function(e, t) {
        const _ = z(n(t)), c = p(_, o.__wbindgen_export, o.__wbindgen_export2), a = g;
        w().setInt32(e + 4, a, !0), w().setInt32(e + 0, c, !0);
      },
      __wbg___wbindgen_is_function_3c846841762788c1: function(e) {
        return typeof n(e) == "function";
      },
      __wbg___wbindgen_is_null_0b605fc6b167c56f: function(e) {
        return n(e) === null;
      },
      __wbg___wbindgen_is_undefined_52709e72fb9f179c: function(e) {
        return n(e) === void 0;
      },
      __wbg___wbindgen_throw_6ddd609b62940d55: function(e, t) {
        throw new Error(u(e, t));
      },
      __wbg__wbg_cb_unref_6b5b6b8576d35cb1: function(e) {
        n(e)._wbg_cb_unref();
      },
      __wbg_beginRenderPass_373f34636d157c43: function() {
        return b(function(e, t) {
          const _ = n(e).beginRenderPass(n(t));
          return i(_);
        }, arguments);
      },
      __wbg_buffer_60b8043cd926067d: function(e) {
        const t = n(e).buffer;
        return i(t);
      },
      __wbg_call_2d781c1f4d5c0ef8: function() {
        return b(function(e, t, _) {
          const c = n(e).call(n(t), n(_));
          return i(c);
        }, arguments);
      },
      __wbg_clientHeight_3d6e452054fdbc3b: function(e) {
        return n(e).clientHeight;
      },
      __wbg_clientWidth_33c7e9c1bcdf4a7e: function(e) {
        return n(e).clientWidth;
      },
      __wbg_configure_b39d6ec9527208fd: function() {
        return b(function(e, t) {
          n(e).configure(n(t));
        }, arguments);
      },
      __wbg_copyTextureToBuffer_f5501895b13306e1: function() {
        return b(function(e, t, _, c) {
          n(e).copyTextureToBuffer(n(t), n(_), n(c));
        }, arguments);
      },
      __wbg_createBindGroupLayout_f5bb5a31b2ac11bf: function() {
        return b(function(e, t) {
          const _ = n(e).createBindGroupLayout(n(t));
          return i(_);
        }, arguments);
      },
      __wbg_createBindGroup_2290306cfa413c74: function(e, t) {
        const _ = n(e).createBindGroup(n(t));
        return i(_);
      },
      __wbg_createBuffer_e2b25dd1471f92f7: function() {
        return b(function(e, t) {
          const _ = n(e).createBuffer(n(t));
          return i(_);
        }, arguments);
      },
      __wbg_createCommandEncoder_80578730e7314357: function(e, t) {
        const _ = n(e).createCommandEncoder(n(t));
        return i(_);
      },
      __wbg_createPipelineLayout_0ef251301bed0c34: function(e, t) {
        const _ = n(e).createPipelineLayout(n(t));
        return i(_);
      },
      __wbg_createRenderPipeline_f9f8aa23f50f8a9c: function() {
        return b(function(e, t) {
          const _ = n(e).createRenderPipeline(n(t));
          return i(_);
        }, arguments);
      },
      __wbg_createSampler_27c37a8245da51a4: function(e, t) {
        const _ = n(e).createSampler(n(t));
        return i(_);
      },
      __wbg_createShaderModule_eb21a131dfb0d4dc: function(e, t) {
        const _ = n(e).createShaderModule(n(t));
        return i(_);
      },
      __wbg_createTexture_284160f981e0075f: function() {
        return b(function(e, t) {
          const _ = n(e).createTexture(n(t));
          return i(_);
        }, arguments);
      },
      __wbg_createView_b09749798973b0f5: function() {
        return b(function(e, t) {
          const _ = n(e).createView(n(t));
          return i(_);
        }, arguments);
      },
      __wbg_debug_4b9b1a2d5972be57: function(e) {
        console.debug(n(e));
      },
      __wbg_document_c0320cd4183c6d9b: function(e) {
        const t = n(e).document;
        return m(t) ? 0 : i(t);
      },
      __wbg_drawIndexed_a60a41b2b0ffdadf: function(e, t, _, c, a, s) {
        n(e).drawIndexed(t >>> 0, _ >>> 0, c >>> 0, a, s >>> 0);
      },
      __wbg_draw_bcc050d6677121b5: function(e, t, _, c, a) {
        n(e).draw(t >>> 0, _ >>> 0, c >>> 0, a >>> 0);
      },
      __wbg_end_c269ebd826210ed1: function(e) {
        n(e).end();
      },
      __wbg_error_8d9a8e04cd1d3588: function(e) {
        console.error(n(e));
      },
      __wbg_error_a6fa202b58aa1cd3: function(e, t) {
        let _, c;
        try {
          _ = e, c = t, console.error(u(e, t));
        } finally {
          o.__wbindgen_export4(_, c, 1);
        }
      },
      __wbg_finish_073e2bc456a4b625: function(e) {
        const t = n(e).finish();
        return i(t);
      },
      __wbg_finish_e43b1b48427f2db0: function(e, t) {
        const _ = n(e).finish(n(t));
        return i(_);
      },
      __wbg_getContext_a9236f98f1f7fe7c: function() {
        return b(function(e, t, _) {
          const c = n(e).getContext(u(t, _));
          return m(c) ? 0 : i(c);
        }, arguments);
      },
      __wbg_getContext_f04bf8f22dcb2d53: function() {
        return b(function(e, t, _) {
          const c = n(e).getContext(u(t, _));
          return m(c) ? 0 : i(c);
        }, arguments);
      },
      __wbg_getCurrentTexture_7edbea16b438c9fc: function() {
        return b(function(e) {
          const t = n(e).getCurrentTexture();
          return i(t);
        }, arguments);
      },
      __wbg_getElementById_d1f25d287b19a833: function(e, t, _) {
        const c = n(e).getElementById(u(t, _));
        return m(c) ? 0 : i(c);
      },
      __wbg_getMappedRange_191c0084744858f0: function() {
        return b(function(e, t, _) {
          const c = n(e).getMappedRange(t, _);
          return i(c);
        }, arguments);
      },
      __wbg_getPreferredCanvasFormat_56e30944cc798353: function(e) {
        const t = n(e).getPreferredCanvasFormat();
        return (y.indexOf(t) + 1 || 96) - 1;
      },
      __wbg_getRandomValues_3f44b700395062e5: function() {
        return b(function(e, t) {
          globalThis.crypto.getRandomValues(P(e, t));
        }, arguments);
      },
      __wbg_get_c7546417fb0bec10: function(e, t) {
        const _ = n(e)[t >>> 0];
        return m(_) ? 0 : i(_);
      },
      __wbg_gpu_7c0927abcc96dd45: function(e) {
        const t = n(e).gpu;
        return i(t);
      },
      __wbg_info_7d4e223bb1a7e671: function(e) {
        console.info(n(e));
      },
      __wbg_instanceof_GpuAdapter_5e451ad6596e2784: function(e) {
        let t;
        try {
          t = n(e) instanceof GPUAdapter;
        } catch {
          t = !1;
        }
        return t;
      },
      __wbg_instanceof_GpuCanvasContext_f70ee27f49f4f884: function(e) {
        let t;
        try {
          t = n(e) instanceof GPUCanvasContext;
        } catch {
          t = !1;
        }
        return t;
      },
      __wbg_instanceof_HtmlCanvasElement_26125339f936be50: function(e) {
        let t;
        try {
          t = n(e) instanceof HTMLCanvasElement;
        } catch {
          t = !1;
        }
        return t;
      },
      __wbg_instanceof_Window_23e677d2c6843922: function(e) {
        let t;
        try {
          t = n(e) instanceof Window;
        } catch {
          t = !1;
        }
        return t;
      },
      __wbg_label_0abc44bf8d3a3e99: function(e, t) {
        const _ = n(t).label, c = p(_, o.__wbindgen_export, o.__wbindgen_export2), a = g;
        w().setInt32(e + 4, a, !0), w().setInt32(e + 0, c, !0);
      },
      __wbg_length_ea16607d7b61445b: function(e) {
        return n(e).length;
      },
      __wbg_limits_764638d29dec49d4: function(e) {
        const t = n(e).limits;
        return i(t);
      },
      __wbg_log_524eedafa26daa59: function(e) {
        console.log(n(e));
      },
      __wbg_mapAsync_1be2f9e8f464f69e: function(e, t, _, c) {
        const a = n(e).mapAsync(t >>> 0, _, c);
        return i(a);
      },
      __wbg_maxBindGroups_c439abd1498fc924: function(e) {
        return n(e).maxBindGroups;
      },
      __wbg_maxBindingsPerBindGroup_186292f383c7b982: function(e) {
        return n(e).maxBindingsPerBindGroup;
      },
      __wbg_maxBufferSize_87b76aa2842d0e8e: function(e) {
        return n(e).maxBufferSize;
      },
      __wbg_maxColorAttachmentBytesPerSample_2ba81ae1e2742413: function(e) {
        return n(e).maxColorAttachmentBytesPerSample;
      },
      __wbg_maxColorAttachments_1ec5191521ef0d22: function(e) {
        return n(e).maxColorAttachments;
      },
      __wbg_maxComputeInvocationsPerWorkgroup_ee67a82206d412d2: function(e) {
        return n(e).maxComputeInvocationsPerWorkgroup;
      },
      __wbg_maxComputeWorkgroupSizeX_0b2b16b802f85a14: function(e) {
        return n(e).maxComputeWorkgroupSizeX;
      },
      __wbg_maxComputeWorkgroupSizeY_00d8aeba9472fdb2: function(e) {
        return n(e).maxComputeWorkgroupSizeY;
      },
      __wbg_maxComputeWorkgroupSizeZ_351fd9dab4c07321: function(e) {
        return n(e).maxComputeWorkgroupSizeZ;
      },
      __wbg_maxComputeWorkgroupStorageSize_881d2b675868eb68: function(e) {
        return n(e).maxComputeWorkgroupStorageSize;
      },
      __wbg_maxComputeWorkgroupsPerDimension_21c223eca6bd6d6b: function(e) {
        return n(e).maxComputeWorkgroupsPerDimension;
      },
      __wbg_maxDynamicStorageBuffersPerPipelineLayout_7155d3f7a514a157: function(e) {
        return n(e).maxDynamicStorageBuffersPerPipelineLayout;
      },
      __wbg_maxDynamicUniformBuffersPerPipelineLayout_76dee9028eaa5322: function(e) {
        return n(e).maxDynamicUniformBuffersPerPipelineLayout;
      },
      __wbg_maxSampledTexturesPerShaderStage_78d018dcd0b999c8: function(e) {
        return n(e).maxSampledTexturesPerShaderStage;
      },
      __wbg_maxSamplersPerShaderStage_0e3ad4d70194a7c2: function(e) {
        return n(e).maxSamplersPerShaderStage;
      },
      __wbg_maxStorageBufferBindingSize_30a1e5c0b8fcd992: function(e) {
        return n(e).maxStorageBufferBindingSize;
      },
      __wbg_maxStorageBuffersPerShaderStage_d77703e9a0d5960e: function(e) {
        return n(e).maxStorageBuffersPerShaderStage;
      },
      __wbg_maxStorageTexturesPerShaderStage_c09e7daf1141067e: function(e) {
        return n(e).maxStorageTexturesPerShaderStage;
      },
      __wbg_maxTextureArrayLayers_44d8badedb4e5245: function(e) {
        return n(e).maxTextureArrayLayers;
      },
      __wbg_maxTextureDimension1D_6d1ff8e56b9cf824: function(e) {
        return n(e).maxTextureDimension1D;
      },
      __wbg_maxTextureDimension2D_5ef5830837d92b7c: function(e) {
        return n(e).maxTextureDimension2D;
      },
      __wbg_maxTextureDimension3D_cfdebbf2b20068cd: function(e) {
        return n(e).maxTextureDimension3D;
      },
      __wbg_maxUniformBufferBindingSize_63dc0c714d2fcebe: function(e) {
        return n(e).maxUniformBufferBindingSize;
      },
      __wbg_maxUniformBuffersPerShaderStage_a52382f8a7dfc816: function(e) {
        return n(e).maxUniformBuffersPerShaderStage;
      },
      __wbg_maxVertexAttributes_4c83ac8c1d442e1c: function(e) {
        return n(e).maxVertexAttributes;
      },
      __wbg_maxVertexBufferArrayStride_955879053ec672f8: function(e) {
        return n(e).maxVertexBufferArrayStride;
      },
      __wbg_maxVertexBuffers_0bb014e62f100c6c: function(e) {
        return n(e).maxVertexBuffers;
      },
      __wbg_minStorageBufferOffsetAlignment_6ed09762e603ac3a: function(e) {
        return n(e).minStorageBufferOffsetAlignment;
      },
      __wbg_minUniformBufferOffsetAlignment_02579f79815cf83c: function(e) {
        return n(e).minUniformBufferOffsetAlignment;
      },
      __wbg_navigator_583ffd4fc14c0f7a: function(e) {
        const t = n(e).navigator;
        return i(t);
      },
      __wbg_navigator_9cebf56f28aa719b: function(e) {
        const t = n(e).navigator;
        return i(t);
      },
      __wbg_new_227d7c05414eb861: function() {
        const e = new Error();
        return i(e);
      },
      __wbg_new_a70fbab9066b301f: function() {
        const e = new Array();
        return i(e);
      },
      __wbg_new_ab79df5bd7c26067: function() {
        const e = new Object();
        return i(e);
      },
      __wbg_new_from_slice_22da9388ac046e50: function(e, t) {
        const _ = new Uint8Array(P(e, t));
        return i(_);
      },
      __wbg_new_typed_aaaeaf29cf802876: function(e, t) {
        try {
          var _ = { a: e, b: t }, c = (s, f) => {
            const l = _.a;
            _.a = 0;
            try {
              return K(l, _.b, s, f);
            } finally {
              _.a = l;
            }
          };
          const a = new Promise(c);
          return i(a);
        } finally {
          _.a = _.b = 0;
        }
      },
      __wbg_new_typed_bccac67128ed885a: function() {
        const e = new Array();
        return i(e);
      },
      __wbg_new_with_byte_offset_and_length_b2ec5bf7b2f35743: function(e, t, _) {
        const c = new Uint8Array(n(e), t >>> 0, _ >>> 0);
        return i(c);
      },
      __wbg_now_e7c6795a7f81e10f: function(e) {
        return n(e).now();
      },
      __wbg_onSubmittedWorkDone_7d532ba1f20a64b3: function(e) {
        const t = n(e).onSubmittedWorkDone();
        return i(t);
      },
      __wbg_performance_3fcf6e32a7e1ed0a: function(e) {
        const t = n(e).performance;
        return i(t);
      },
      __wbg_prototypesetcall_d62e5099504357e6: function(e, t, _) {
        Uint8Array.prototype.set.call(P(e, t), n(_));
      },
      __wbg_push_e87b0e732085a946: function(e, t) {
        return n(e).push(n(t));
      },
      __wbg_querySelectorAll_ccbf0696a1c6fed8: function() {
        return b(function(e, t, _) {
          const c = n(e).querySelectorAll(u(t, _));
          return i(c);
        }, arguments);
      },
      __wbg_queueMicrotask_0c399741342fb10f: function(e) {
        const t = n(e).queueMicrotask;
        return i(t);
      },
      __wbg_queueMicrotask_a082d78ce798393e: function(e) {
        queueMicrotask(n(e));
      },
      __wbg_queue_5eda23116e5d3adb: function(e) {
        const t = n(e).queue;
        return i(t);
      },
      __wbg_requestAdapter_8efca1b953fd13aa: function(e, t) {
        const _ = n(e).requestAdapter(n(t));
        return i(_);
      },
      __wbg_requestDevice_290c73161fe959d5: function(e, t) {
        const _ = n(e).requestDevice(n(t));
        return i(_);
      },
      __wbg_resolve_ae8d83246e5bcc12: function(e) {
        const t = Promise.resolve(n(e));
        return i(t);
      },
      __wbg_setBindGroup_29f4a44dff76f1a4: function(e, t, _) {
        n(e).setBindGroup(t >>> 0, n(_));
      },
      __wbg_setBindGroup_35a4830ac2c27742: function() {
        return b(function(e, t, _, c, a, s, f) {
          n(e).setBindGroup(t >>> 0, n(_), we(c, a), s, f >>> 0);
        }, arguments);
      },
      __wbg_setIndexBuffer_924197dc97dbb679: function(e, t, _, c, a) {
        n(e).setIndexBuffer(n(t), G[_], c, a);
      },
      __wbg_setIndexBuffer_a400322dea5437f7: function(e, t, _, c) {
        n(e).setIndexBuffer(n(t), G[_], c);
      },
      __wbg_setPipeline_e6ea6756d71b19a7: function(e, t) {
        n(e).setPipeline(n(t));
      },
      __wbg_setScissorRect_eeb4f61d4b860d7a: function(e, t, _, c, a) {
        n(e).setScissorRect(t >>> 0, _ >>> 0, c >>> 0, a >>> 0);
      },
      __wbg_setVertexBuffer_58f30a4873b36907: function(e, t, _, c) {
        n(e).setVertexBuffer(t >>> 0, n(_), c);
      },
      __wbg_setVertexBuffer_7aa508f017477005: function(e, t, _, c, a) {
        n(e).setVertexBuffer(t >>> 0, n(_), c, a);
      },
      __wbg_setViewport_014b4c4d1101ba6b: function(e, t, _, c, a, s, f) {
        n(e).setViewport(t, _, c, a, s, f);
      },
      __wbg_set_282384002438957f: function(e, t, _) {
        n(e)[t >>> 0] = d(_);
      },
      __wbg_set_6be42768c690e380: function(e, t, _) {
        n(e)[d(t)] = d(_);
      },
      __wbg_set_7eaa4f96924fd6b3: function() {
        return b(function(e, t, _) {
          return Reflect.set(n(e), n(t), n(_));
        }, arguments);
      },
      __wbg_set_a_6f1653ca7319cdcf: function(e, t) {
        n(e).a = t;
      },
      __wbg_set_access_cbee993a36feed10: function(e, t) {
        n(e).access = ae[t];
      },
      __wbg_set_address_mode_u_38e255cd89ce1977: function(e, t) {
        n(e).addressModeU = O[t];
      },
      __wbg_set_address_mode_v_513f843d6e3c9dbd: function(e, t) {
        n(e).addressModeV = O[t];
      },
      __wbg_set_address_mode_w_801f70901a90ed5a: function(e, t) {
        n(e).addressModeW = O[t];
      },
      __wbg_set_alpha_0a28ffc800461787: function(e, t) {
        n(e).alpha = n(t);
      },
      __wbg_set_alpha_mode_55b4f33e93691fe8: function(e, t) {
        n(e).alphaMode = te[t];
      },
      __wbg_set_alpha_to_coverage_enabled_ec44695cc0d0e961: function(e, t) {
        n(e).alphaToCoverageEnabled = t !== 0;
      },
      __wbg_set_array_layer_count_e774b6d4a5334e63: function(e, t) {
        n(e).arrayLayerCount = t >>> 0;
      },
      __wbg_set_array_stride_11c840b41b728354: function(e, t) {
        n(e).arrayStride = t;
      },
      __wbg_set_aspect_2503cdfcdcc17373: function(e, t) {
        n(e).aspect = se[t];
      },
      __wbg_set_attributes_ac1030b589bf253a: function(e, t) {
        n(e).attributes = n(t);
      },
      __wbg_set_b_d5b23064b0492744: function(e, t) {
        n(e).b = t;
      },
      __wbg_set_base_array_layer_f64cdadf250d1a9b: function(e, t) {
        n(e).baseArrayLayer = t >>> 0;
      },
      __wbg_set_base_mip_level_74fc97c2aaf8fc33: function(e, t) {
        n(e).baseMipLevel = t >>> 0;
      },
      __wbg_set_beginning_of_pass_write_index_348e7f2f53a86db0: function(e, t) {
        n(e).beginningOfPassWriteIndex = t >>> 0;
      },
      __wbg_set_bind_group_layouts_6f13eb021a550053: function(e, t) {
        n(e).bindGroupLayouts = n(t);
      },
      __wbg_set_binding_2240d98479c0c256: function(e, t) {
        n(e).binding = t >>> 0;
      },
      __wbg_set_binding_5296904f2a4c7e25: function(e, t) {
        n(e).binding = t >>> 0;
      },
      __wbg_set_blend_4aea897cd7d3c0f8: function(e, t) {
        n(e).blend = n(t);
      },
      __wbg_set_buffer_2e7d1f7814caf92b: function(e, t) {
        n(e).buffer = n(t);
      },
      __wbg_set_buffer_ba8ed06078d347f7: function(e, t) {
        n(e).buffer = n(t);
      },
      __wbg_set_buffer_fc9285180932669f: function(e, t) {
        n(e).buffer = n(t);
      },
      __wbg_set_buffers_72754529595d4bc0: function(e, t) {
        n(e).buffers = n(t);
      },
      __wbg_set_bytes_per_row_9425e8d6a11b52dc: function(e, t) {
        n(e).bytesPerRow = t >>> 0;
      },
      __wbg_set_clear_value_1171de96edbc21fe: function(e, t) {
        n(e).clearValue = n(t);
      },
      __wbg_set_code_27a25a855d3fbc6d: function(e, t, _) {
        n(e).code = u(t, _);
      },
      __wbg_set_color_attachments_4516b6dfb4ad987b: function(e, t) {
        n(e).colorAttachments = n(t);
      },
      __wbg_set_color_f2ac28bdc576c010: function(e, t) {
        n(e).color = n(t);
      },
      __wbg_set_compare_2c8ee8ccaa2b6b5d: function(e, t) {
        n(e).compare = W[t];
      },
      __wbg_set_compare_cbf49b43d3211833: function(e, t) {
        n(e).compare = W[t];
      },
      __wbg_set_count_53854513da5c0e04: function(e, t) {
        n(e).count = t >>> 0;
      },
      __wbg_set_cull_mode_3852dd4cff56dd90: function(e, t) {
        n(e).cullMode = ne[t];
      },
      __wbg_set_depth_bias_c20861a58fc2b8d9: function(e, t) {
        n(e).depthBias = t;
      },
      __wbg_set_depth_bias_clamp_eecc04d702f9402e: function(e, t) {
        n(e).depthBiasClamp = t;
      },
      __wbg_set_depth_bias_slope_scale_b2a251d3d4c65018: function(e, t) {
        n(e).depthBiasSlopeScale = t;
      },
      __wbg_set_depth_clear_value_fca9e379a0cdff8f: function(e, t) {
        n(e).depthClearValue = t;
      },
      __wbg_set_depth_compare_7883e52aad39b925: function(e, t) {
        n(e).depthCompare = W[t];
      },
      __wbg_set_depth_fail_op_1d11c8e03d061484: function(e, t) {
        n(e).depthFailOp = D[t];
      },
      __wbg_set_depth_load_op_7e95e67c69e09c5e: function(e, t) {
        n(e).depthLoadOp = M[t];
      },
      __wbg_set_depth_or_array_layers_36ef1df107b6b651: function(e, t) {
        n(e).depthOrArrayLayers = t >>> 0;
      },
      __wbg_set_depth_read_only_0c5e726b56520b08: function(e, t) {
        n(e).depthReadOnly = t !== 0;
      },
      __wbg_set_depth_stencil_17e2d1710f4e07ae: function(e, t) {
        n(e).depthStencil = n(t);
      },
      __wbg_set_depth_stencil_attachment_a7b5eca74b7ddcfb: function(e, t) {
        n(e).depthStencilAttachment = n(t);
      },
      __wbg_set_depth_store_op_1b4cc257f121a4e7: function(e, t) {
        n(e).depthStoreOp = I[t];
      },
      __wbg_set_depth_write_enabled_1551f99ae66d959e: function(e, t) {
        n(e).depthWriteEnabled = t !== 0;
      },
      __wbg_set_device_846227515bb0301a: function(e, t) {
        n(e).device = n(t);
      },
      __wbg_set_dimension_7454baa9c745cf06: function(e, t) {
        n(e).dimension = ue[t];
      },
      __wbg_set_dimension_9d314669636abc65: function(e, t) {
        n(e).dimension = F[t];
      },
      __wbg_set_dst_factor_8397030245674624: function(e, t) {
        n(e).dstFactor = R[t];
      },
      __wbg_set_e80615d7a9a43981: function(e, t, _) {
        n(e).set(n(t), _ >>> 0);
      },
      __wbg_set_end_of_pass_write_index_4600a261d0317ecb: function(e, t) {
        n(e).endOfPassWriteIndex = t >>> 0;
      },
      __wbg_set_entries_4d13c932343146c3: function(e, t) {
        n(e).entries = n(t);
      },
      __wbg_set_entries_7e6b569918b11bf4: function(e, t) {
        n(e).entries = n(t);
      },
      __wbg_set_entry_point_7248ed25fb9070c7: function(e, t, _) {
        n(e).entryPoint = u(t, _);
      },
      __wbg_set_entry_point_b01eb3970a1dcb95: function(e, t, _) {
        n(e).entryPoint = u(t, _);
      },
      __wbg_set_external_texture_cf6cf39036321145: function(e, t) {
        n(e).externalTexture = n(t);
      },
      __wbg_set_fail_op_ac8f2b4c077715b1: function(e, t) {
        n(e).failOp = D[t];
      },
      __wbg_set_format_12bcbdd3428cd4b5: function(e, t) {
        n(e).format = y[t];
      },
      __wbg_set_format_1fc8a436841b29c8: function(e, t) {
        n(e).format = y[t];
      },
      __wbg_set_format_2a42ed14de233ae5: function(e, t) {
        n(e).format = be[t];
      },
      __wbg_set_format_3759d043ddc658d4: function(e, t) {
        n(e).format = y[t];
      },
      __wbg_set_format_b08e529cc1612d7b: function(e, t) {
        n(e).format = y[t];
      },
      __wbg_set_format_e0cf5a237864edb6: function(e, t) {
        n(e).format = y[t];
      },
      __wbg_set_format_ffa0a97f114a945a: function(e, t) {
        n(e).format = y[t];
      },
      __wbg_set_fragment_703ddd6f5db6e4af: function(e, t) {
        n(e).fragment = n(t);
      },
      __wbg_set_front_face_17a3723085696d9a: function(e, t) {
        n(e).frontFace = _e[t];
      },
      __wbg_set_g_4cc3b3e3231ca6f8: function(e, t) {
        n(e).g = t;
      },
      __wbg_set_has_dynamic_offset_dc25aba64b9bd3ff: function(e, t) {
        n(e).hasDynamicOffset = t !== 0;
      },
      __wbg_set_height_98a1a397672657e2: function(e, t) {
        n(e).height = t >>> 0;
      },
      __wbg_set_height_ac705ece3aa08c95: function(e, t) {
        n(e).height = t >>> 0;
      },
      __wbg_set_height_b6548a01bdcb689a: function(e, t) {
        n(e).height = t >>> 0;
      },
      __wbg_set_label_10bd19b972ff1ba6: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_16cff4ff3c381368: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_343ceab4761679d7: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_403725ced930414e: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_62b82f9361718fb9: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_7d448e8a777d0d37: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_900e563567315063: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_98bef61fcbcecdde: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_9d2ce197e447a967: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_b5d7ff5f8e4fbaac: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_ba288fbac1259847: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_e1bd2437f39d21f3: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_label_e4debe6dc9ea319b: function(e, t, _) {
        n(e).label = u(t, _);
      },
      __wbg_set_layout_53be3643dc5dbbbe: function(e, t) {
        n(e).layout = n(t);
      },
      __wbg_set_layout_ca5f863d331bb6b4: function(e, t) {
        n(e).layout = n(t);
      },
      __wbg_set_load_op_91d2cbf2912c96fd: function(e, t) {
        n(e).loadOp = M[t];
      },
      __wbg_set_lod_max_clamp_01800ff5df00cc8e: function(e, t) {
        n(e).lodMaxClamp = t;
      },
      __wbg_set_lod_min_clamp_fe71be084b04bd97: function(e, t) {
        n(e).lodMinClamp = t;
      },
      __wbg_set_mag_filter_a6df09d1943d5caa: function(e, t) {
        n(e).magFilter = V[t];
      },
      __wbg_set_mapped_at_creation_eb954cf5fdb9bc25: function(e, t) {
        n(e).mappedAtCreation = t !== 0;
      },
      __wbg_set_mask_47a41aae6631771f: function(e, t) {
        n(e).mask = t >>> 0;
      },
      __wbg_set_max_anisotropy_418bd200a56097a0: function(e, t) {
        n(e).maxAnisotropy = t;
      },
      __wbg_set_min_binding_size_d0315b751370234c: function(e, t) {
        n(e).minBindingSize = t;
      },
      __wbg_set_min_filter_5b27a7eb3f5ea88a: function(e, t) {
        n(e).minFilter = V[t];
      },
      __wbg_set_mip_level_b50dccbd04935c98: function(e, t) {
        n(e).mipLevel = t >>> 0;
      },
      __wbg_set_mip_level_count_307eb64d9d29e3a6: function(e, t) {
        n(e).mipLevelCount = t >>> 0;
      },
      __wbg_set_mip_level_count_fe7f73daa6021aaa: function(e, t) {
        n(e).mipLevelCount = t >>> 0;
      },
      __wbg_set_mipmap_filter_e1543204e8199db0: function(e, t) {
        n(e).mipmapFilter = re[t];
      },
      __wbg_set_module_9afd1b80ff72cee9: function(e, t) {
        n(e).module = n(t);
      },
      __wbg_set_module_ffe8f8e909e9fdcf: function(e, t) {
        n(e).module = n(t);
      },
      __wbg_set_multisample_957afdd96685c6f5: function(e, t) {
        n(e).multisample = n(t);
      },
      __wbg_set_multisampled_84e304d3a68838ea: function(e, t) {
        n(e).multisampled = t !== 0;
      },
      __wbg_set_offset_157c6bc4fd6ec4b1: function(e, t) {
        n(e).offset = t;
      },
      __wbg_set_offset_3e78f3e530cf8049: function(e, t) {
        n(e).offset = t;
      },
      __wbg_set_offset_bea112c360dc7f2b: function(e, t) {
        n(e).offset = t;
      },
      __wbg_set_operation_6c5fd88df90bc7b2: function(e, t) {
        n(e).operation = Q[t];
      },
      __wbg_set_origin_dec4f4c36f9f79f6: function(e, t) {
        n(e).origin = n(t);
      },
      __wbg_set_pass_op_461dabd5ee4ea1b7: function(e, t) {
        n(e).passOp = D[t];
      },
      __wbg_set_power_preference_a4ce891b22ea2b05: function(e, t) {
        n(e).powerPreference = ce[t];
      },
      __wbg_set_primitive_eb8abbc5e7f278a4: function(e, t) {
        n(e).primitive = n(t);
      },
      __wbg_set_query_set_849fb32875f137d7: function(e, t) {
        n(e).querySet = n(t);
      },
      __wbg_set_r_5fa0f548248c394c: function(e, t) {
        n(e).r = t;
      },
      __wbg_set_required_features_98a83c7003fd73d5: function(e, t) {
        n(e).requiredFeatures = n(t);
      },
      __wbg_set_resolve_target_1ff405e060e2d32e: function(e, t) {
        n(e).resolveTarget = n(t);
      },
      __wbg_set_resource_1409c14d4d6b5a50: function(e, t) {
        n(e).resource = n(t);
      },
      __wbg_set_rows_per_image_9cfda8920e669db0: function(e, t) {
        n(e).rowsPerImage = t >>> 0;
      },
      __wbg_set_sample_count_95a9892a60894677: function(e, t) {
        n(e).sampleCount = t >>> 0;
      },
      __wbg_set_sample_type_f8f7b39d62e7b29c: function(e, t) {
        n(e).sampleType = fe[t];
      },
      __wbg_set_sampler_a2277e90dfe7395f: function(e, t) {
        n(e).sampler = n(t);
      },
      __wbg_set_shader_location_cdbcf5cf84a6cbcb: function(e, t) {
        n(e).shaderLocation = t >>> 0;
      },
      __wbg_set_size_6f271c4c28c18e1b: function(e, t) {
        n(e).size = n(t);
      },
      __wbg_set_size_7ec162511b3bad1f: function(e, t) {
        n(e).size = t;
      },
      __wbg_set_size_ca765d983baccefd: function(e, t) {
        n(e).size = t;
      },
      __wbg_set_src_factor_e96f05a25f8383ed: function(e, t) {
        n(e).srcFactor = R[t];
      },
      __wbg_set_stencil_back_5c8971274cbcddcf: function(e, t) {
        n(e).stencilBack = n(t);
      },
      __wbg_set_stencil_clear_value_89ba97b367fa1385: function(e, t) {
        n(e).stencilClearValue = t >>> 0;
      },
      __wbg_set_stencil_front_69f85bf4a6f02cb2: function(e, t) {
        n(e).stencilFront = n(t);
      },
      __wbg_set_stencil_load_op_a3e2c3a6f20d4da5: function(e, t) {
        n(e).stencilLoadOp = M[t];
      },
      __wbg_set_stencil_read_mask_86a08afb2665c29b: function(e, t) {
        n(e).stencilReadMask = t >>> 0;
      },
      __wbg_set_stencil_read_only_dd058fe8c6a1f6ae: function(e, t) {
        n(e).stencilReadOnly = t !== 0;
      },
      __wbg_set_stencil_store_op_87c97415636844c9: function(e, t) {
        n(e).stencilStoreOp = I[t];
      },
      __wbg_set_stencil_write_mask_7844d8a057a87a58: function(e, t) {
        n(e).stencilWriteMask = t >>> 0;
      },
      __wbg_set_step_mode_285f2e428148f3b4: function(e, t) {
        n(e).stepMode = de[t];
      },
      __wbg_set_storage_texture_373b9fc0e534dd33: function(e, t) {
        n(e).storageTexture = n(t);
      },
      __wbg_set_store_op_94575f47253d270d: function(e, t) {
        n(e).storeOp = I[t];
      },
      __wbg_set_strip_index_format_aeb7aa0e95e6285d: function(e, t) {
        n(e).stripIndexFormat = G[t];
      },
      __wbg_set_targets_93553735385af349: function(e, t) {
        n(e).targets = n(t);
      },
      __wbg_set_texture_6003a9e79918bf8a: function(e, t) {
        n(e).texture = n(t);
      },
      __wbg_set_texture_c5a457625c071b25: function(e, t) {
        n(e).texture = n(t);
      },
      __wbg_set_timestamp_writes_0603b32a31ee6205: function(e, t) {
        n(e).timestampWrites = n(t);
      },
      __wbg_set_topology_5e4eb809635ea291: function(e, t) {
        n(e).topology = oe[t];
      },
      __wbg_set_type_0e707d4c06fc2b7b: function(e, t) {
        n(e).type = ie[t];
      },
      __wbg_set_type_6fe4c5f460401ee0: function(e, t) {
        n(e).type = ee[t];
      },
      __wbg_set_unclipped_depth_e9a2451e4fa0277a: function(e, t) {
        n(e).unclippedDepth = t !== 0;
      },
      __wbg_set_usage_5abd566becc087bb: function(e, t) {
        n(e).usage = t >>> 0;
      },
      __wbg_set_usage_61967f166fba5e13: function(e, t) {
        n(e).usage = t >>> 0;
      },
      __wbg_set_usage_d0a75d4429098a06: function(e, t) {
        n(e).usage = t >>> 0;
      },
      __wbg_set_usage_f0bb325677668e77: function(e, t) {
        n(e).usage = t >>> 0;
      },
      __wbg_set_vertex_2525cfcd959b2add: function(e, t) {
        n(e).vertex = n(t);
      },
      __wbg_set_view_57d232eea19739c3: function(e, t) {
        n(e).view = n(t);
      },
      __wbg_set_view_dimension_49cfda500f1dea55: function(e, t) {
        n(e).viewDimension = F[t];
      },
      __wbg_set_view_dimension_a669c29ec3b0813a: function(e, t) {
        n(e).viewDimension = F[t];
      },
      __wbg_set_view_ffadd767d5e9b839: function(e, t) {
        n(e).view = n(t);
      },
      __wbg_set_view_formats_70a1fcabcd34282a: function(e, t) {
        n(e).viewFormats = n(t);
      },
      __wbg_set_view_formats_83865b9cdfda5cb6: function(e, t) {
        n(e).viewFormats = n(t);
      },
      __wbg_set_visibility_088046ee77c33b1d: function(e, t) {
        n(e).visibility = t >>> 0;
      },
      __wbg_set_width_576343a4a7f2cf28: function(e, t) {
        n(e).width = t >>> 0;
      },
      __wbg_set_width_c0fcaa2da53cd540: function(e, t) {
        n(e).width = t >>> 0;
      },
      __wbg_set_width_e96e07f8255ad913: function(e, t) {
        n(e).width = t >>> 0;
      },
      __wbg_set_write_mask_76041c03688571cd: function(e, t) {
        n(e).writeMask = t >>> 0;
      },
      __wbg_set_x_fdd6aca9a2390926: function(e, t) {
        n(e).x = t >>> 0;
      },
      __wbg_set_y_410a18c5811abf4c: function(e, t) {
        n(e).y = t >>> 0;
      },
      __wbg_set_z_f7f1ae8afd3a9308: function(e, t) {
        n(e).z = t >>> 0;
      },
      __wbg_stack_3b0d974bbf31e44f: function(e, t) {
        const _ = n(t).stack, c = p(_, o.__wbindgen_export, o.__wbindgen_export2), a = g;
        w().setInt32(e + 4, a, !0), w().setInt32(e + 0, c, !0);
      },
      __wbg_static_accessor_GLOBAL_8adb955bd33fac2f: function() {
        const e = typeof global > "u" ? null : global;
        return m(e) ? 0 : i(e);
      },
      __wbg_static_accessor_GLOBAL_THIS_ad356e0db91c7913: function() {
        const e = typeof globalThis > "u" ? null : globalThis;
        return m(e) ? 0 : i(e);
      },
      __wbg_static_accessor_SELF_f207c857566db248: function() {
        const e = typeof self > "u" ? null : self;
        return m(e) ? 0 : i(e);
      },
      __wbg_static_accessor_WINDOW_bb9f1ba69d61b386: function() {
        const e = typeof window > "u" ? null : window;
        return m(e) ? 0 : i(e);
      },
      __wbg_submit_21302eebe551e30d: function(e, t) {
        n(e).submit(n(t));
      },
      __wbg_then_098abe61755d12f6: function(e, t) {
        const _ = n(e).then(n(t));
        return i(_);
      },
      __wbg_then_9e335f6dd892bc11: function(e, t, _) {
        const c = n(e).then(n(t), n(_));
        return i(c);
      },
      __wbg_then_bc59d1943397ca4e: function(e, t, _) {
        const c = n(e).then(n(t), n(_));
        return i(c);
      },
      __wbg_unmap_b819b8b402db13cc: function(e) {
        n(e).unmap();
      },
      __wbg_warn_69424c2d92a2fa73: function(e) {
        console.warn(n(e));
      },
      __wbg_webviewer_new: function(e) {
        const t = T.__wrap(e);
        return i(t);
      },
      __wbg_writeBuffer_c6919ed0c4aaeef5: function() {
        return b(function(e, t, _, c, a, s) {
          n(e).writeBuffer(n(t), _, n(c), a, s);
        }, arguments);
      },
      __wbindgen_cast_0000000000000001: function(e, t) {
        const _ = j(e, t, o.__wasm_bindgen_func_elem_9952, Z);
        return i(_);
      },
      __wbindgen_cast_0000000000000002: function(e, t) {
        const _ = j(e, t, o.__wasm_bindgen_func_elem_10886, J);
        return i(_);
      },
      __wbindgen_cast_0000000000000003: function(e) {
        return i(e);
      },
      __wbindgen_cast_0000000000000004: function(e, t) {
        const _ = P(e, t);
        return i(_);
      },
      __wbindgen_cast_0000000000000005: function(e, t) {
        const _ = u(e, t);
        return i(_);
      },
      __wbindgen_cast_0000000000000006: function(e) {
        const t = BigInt.asUintN(64, e);
        return i(t);
      },
      __wbindgen_object_clone_ref: function(e) {
        const t = n(e);
        return i(t);
      },
      __wbindgen_object_drop_ref: function(e) {
        d(e);
      }
    }
  };
}
function Z(r, e, t) {
  o.__wasm_bindgen_func_elem_10007(r, e, i(t));
}
function J(r, e, t) {
  try {
    const a = o.__wbindgen_add_to_stack_pointer(-16);
    o.__wasm_bindgen_func_elem_11472(a, r, e, i(t));
    var _ = w().getInt32(a + 0, !0), c = w().getInt32(a + 4, !0);
    if (c)
      throw d(_);
  } finally {
    o.__wbindgen_add_to_stack_pointer(16);
  }
}
function K(r, e, t, _) {
  o.__wasm_bindgen_func_elem_11485(r, e, i(t), i(_));
}
const O = ["clamp-to-edge", "repeat", "mirror-repeat"], R = ["zero", "one", "src", "one-minus-src", "src-alpha", "one-minus-src-alpha", "dst", "one-minus-dst", "dst-alpha", "one-minus-dst-alpha", "src-alpha-saturated", "constant", "one-minus-constant", "src1", "one-minus-src1", "src1-alpha", "one-minus-src1-alpha"], Q = ["add", "subtract", "reverse-subtract", "min", "max"], ee = ["uniform", "storage", "read-only-storage"], te = ["opaque", "premultiplied"], W = ["never", "less", "equal", "less-equal", "greater", "not-equal", "greater-equal", "always"], ne = ["none", "front", "back"], V = ["nearest", "linear"], _e = ["ccw", "cw"], G = ["uint16", "uint32"], M = ["load", "clear"], re = ["nearest", "linear"], ce = ["low-power", "high-performance"], oe = ["point-list", "line-list", "line-strip", "triangle-list", "triangle-strip"], ie = ["filtering", "non-filtering", "comparison"], D = ["keep", "zero", "replace", "invert", "increment-clamp", "decrement-clamp", "increment-wrap", "decrement-wrap"], ae = ["write-only", "read-only", "read-write"], I = ["store", "discard"], se = ["all", "stencil-only", "depth-only"], ue = ["1d", "2d", "3d"], y = ["r8unorm", "r8snorm", "r8uint", "r8sint", "r16uint", "r16sint", "r16float", "rg8unorm", "rg8snorm", "rg8uint", "rg8sint", "r32uint", "r32sint", "r32float", "rg16uint", "rg16sint", "rg16float", "rgba8unorm", "rgba8unorm-srgb", "rgba8snorm", "rgba8uint", "rgba8sint", "bgra8unorm", "bgra8unorm-srgb", "rgb9e5ufloat", "rgb10a2uint", "rgb10a2unorm", "rg11b10ufloat", "rg32uint", "rg32sint", "rg32float", "rgba16uint", "rgba16sint", "rgba16float", "rgba32uint", "rgba32sint", "rgba32float", "stencil8", "depth16unorm", "depth24plus", "depth24plus-stencil8", "depth32float", "depth32float-stencil8", "bc1-rgba-unorm", "bc1-rgba-unorm-srgb", "bc2-rgba-unorm", "bc2-rgba-unorm-srgb", "bc3-rgba-unorm", "bc3-rgba-unorm-srgb", "bc4-r-unorm", "bc4-r-snorm", "bc5-rg-unorm", "bc5-rg-snorm", "bc6h-rgb-ufloat", "bc6h-rgb-float", "bc7-rgba-unorm", "bc7-rgba-unorm-srgb", "etc2-rgb8unorm", "etc2-rgb8unorm-srgb", "etc2-rgb8a1unorm", "etc2-rgb8a1unorm-srgb", "etc2-rgba8unorm", "etc2-rgba8unorm-srgb", "eac-r11unorm", "eac-r11snorm", "eac-rg11unorm", "eac-rg11snorm", "astc-4x4-unorm", "astc-4x4-unorm-srgb", "astc-5x4-unorm", "astc-5x4-unorm-srgb", "astc-5x5-unorm", "astc-5x5-unorm-srgb", "astc-6x5-unorm", "astc-6x5-unorm-srgb", "astc-6x6-unorm", "astc-6x6-unorm-srgb", "astc-8x5-unorm", "astc-8x5-unorm-srgb", "astc-8x6-unorm", "astc-8x6-unorm-srgb", "astc-8x8-unorm", "astc-8x8-unorm-srgb", "astc-10x5-unorm", "astc-10x5-unorm-srgb", "astc-10x6-unorm", "astc-10x6-unorm-srgb", "astc-10x8-unorm", "astc-10x8-unorm-srgb", "astc-10x10-unorm", "astc-10x10-unorm-srgb", "astc-12x10-unorm", "astc-12x10-unorm-srgb", "astc-12x12-unorm", "astc-12x12-unorm-srgb"], fe = ["float", "unfilterable-float", "depth", "sint", "uint"], F = ["1d", "2d", "2d-array", "cube", "cube-array", "3d"], be = ["uint8", "uint8x2", "uint8x4", "sint8", "sint8x2", "sint8x4", "unorm8", "unorm8x2", "unorm8x4", "snorm8", "snorm8x2", "snorm8x4", "uint16", "uint16x2", "uint16x4", "sint16", "sint16x2", "sint16x4", "unorm16", "unorm16x2", "unorm16x4", "snorm16", "snorm16x2", "snorm16x4", "float16", "float16x2", "float16x4", "float32", "float32x2", "float32x3", "float32x4", "uint32", "uint32x2", "uint32x3", "uint32x4", "sint32", "sint32x2", "sint32x3", "sint32x4", "unorm10-10-10-2", "unorm8x4-bgra"], de = ["vertex", "instance"], U = typeof FinalizationRegistry > "u" ? { register: () => {
}, unregister: () => {
} } : new FinalizationRegistry((r) => o.__wbg_webviewer_free(r >>> 0, 1));
function i(r) {
  A === x.length && x.push(x.length + 1);
  const e = A;
  return A = x[e], x[e] = r, e;
}
const q = typeof FinalizationRegistry > "u" ? { register: () => {
}, unregister: () => {
} } : new FinalizationRegistry((r) => r.dtor(r.a, r.b));
function z(r) {
  const e = typeof r;
  if (e == "number" || e == "boolean" || r == null)
    return `${r}`;
  if (e == "string")
    return `"${r}"`;
  if (e == "symbol") {
    const c = r.description;
    return c == null ? "Symbol" : `Symbol(${c})`;
  }
  if (e == "function") {
    const c = r.name;
    return typeof c == "string" && c.length > 0 ? `Function(${c})` : "Function";
  }
  if (Array.isArray(r)) {
    const c = r.length;
    let a = "[";
    c > 0 && (a += z(r[0]));
    for (let s = 1; s < c; s++)
      a += ", " + z(r[s]);
    return a += "]", a;
  }
  const t = /\[object ([^\]]+)\]/.exec(toString.call(r));
  let _;
  if (t && t.length > 1)
    _ = t[1];
  else
    return toString.call(r);
  if (_ == "Object")
    try {
      return "Object(" + JSON.stringify(r) + ")";
    } catch {
      return "Object";
    }
  return r instanceof Error ? `${r.name}: ${r.message}
${r.stack}` : _;
}
function ge(r) {
  r < 1028 || (x[r] = A, A = r);
}
function we(r, e) {
  return r = r >>> 0, le().subarray(r / 4, r / 4 + e);
}
function P(r, e) {
  return r = r >>> 0, S().subarray(r / 1, r / 1 + e);
}
let h = null;
function w() {
  return (h === null || h.buffer.detached === !0 || h.buffer.detached === void 0 && h.buffer !== o.memory.buffer) && (h = new DataView(o.memory.buffer)), h;
}
function u(r, e) {
  return r = r >>> 0, xe(r, e);
}
let v = null;
function le() {
  return (v === null || v.byteLength === 0) && (v = new Uint32Array(o.memory.buffer)), v;
}
let B = null;
function S() {
  return (B === null || B.byteLength === 0) && (B = new Uint8Array(o.memory.buffer)), B;
}
function n(r) {
  return x[r];
}
function b(r, e) {
  try {
    return r.apply(this, e);
  } catch (t) {
    o.__wbindgen_export3(i(t));
  }
}
let x = new Array(1024).fill(void 0);
x.push(void 0, null, !0, !1);
let A = x.length;
function m(r) {
  return r == null;
}
function j(r, e, t, _) {
  const c = { a: r, b: e, cnt: 1, dtor: t }, a = (...s) => {
    c.cnt++;
    const f = c.a;
    c.a = 0;
    try {
      return _(f, c.b, ...s);
    } finally {
      c.a = f, a._wbg_cb_unref();
    }
  };
  return a._wbg_cb_unref = () => {
    --c.cnt === 0 && (c.dtor(c.a, c.b), c.a = 0, q.unregister(c));
  }, q.register(a, c, c), a;
}
function me(r, e) {
  const t = e(r.length * 1, 1) >>> 0;
  return S().set(r, t / 1), g = r.length, t;
}
function p(r, e, t) {
  if (t === void 0) {
    const f = C.encode(r), l = e(f.length, 1) >>> 0;
    return S().subarray(l, l + f.length).set(f), g = f.length, l;
  }
  let _ = r.length, c = e(_, 1) >>> 0;
  const a = S();
  let s = 0;
  for (; s < _; s++) {
    const f = r.charCodeAt(s);
    if (f > 127) break;
    a[c + s] = f;
  }
  if (s !== _) {
    s !== 0 && (r = r.slice(s)), c = t(c, _, _ = s + r.length * 3, 1) >>> 0;
    const f = S().subarray(c + s, c + _), l = C.encodeInto(r, f);
    s += l.written, c = t(c, _, s, 1) >>> 0;
  }
  return g = s, c;
}
function d(r) {
  const e = n(r);
  return ge(r), e;
}
let k = new TextDecoder("utf-8", { ignoreBOM: !0, fatal: !0 });
k.decode();
const pe = 2146435072;
let L = 0;
function xe(r, e) {
  return L += e, L >= pe && (k = new TextDecoder("utf-8", { ignoreBOM: !0, fatal: !0 }), k.decode(), L = e), k.decode(S().subarray(r, r + e));
}
const C = new TextEncoder();
"encodeInto" in C || (C.encodeInto = function(r, e) {
  const t = C.encode(r);
  return e.set(t), {
    read: r.length,
    written: t.length
  };
});
let g = 0, o;
function $(r, e) {
  return o = r.exports, h = null, v = null, B = null, o.__wbindgen_start(), o;
}
async function ye(r, e) {
  if (typeof Response == "function" && r instanceof Response) {
    if (typeof WebAssembly.instantiateStreaming == "function")
      try {
        return await WebAssembly.instantiateStreaming(r, e);
      } catch (c) {
        if (r.ok && t(r.type) && r.headers.get("Content-Type") !== "application/wasm")
          console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", c);
        else
          throw c;
      }
    const _ = await r.arrayBuffer();
    return await WebAssembly.instantiate(_, e);
  } else {
    const _ = await WebAssembly.instantiate(r, e);
    return _ instanceof WebAssembly.Instance ? { instance: _, module: r } : _;
  }
  function t(_) {
    switch (_) {
      case "basic":
      case "cors":
      case "default":
        return !0;
    }
    return !1;
  }
}
function Se(r) {
  if (o !== void 0) return o;
  r !== void 0 && (Object.getPrototypeOf(r) === Object.prototype ? { module: r } = r : console.warn("using deprecated parameters for `initSync()`; pass a single object instead"));
  const e = E();
  r instanceof WebAssembly.Module || (r = new WebAssembly.Module(r));
  const t = new WebAssembly.Instance(r, e);
  return $(t);
}
async function ve(r) {
  if (o !== void 0) return o;
  r !== void 0 && (Object.getPrototypeOf(r) === Object.prototype ? { module_or_path: r } = r : console.warn("using deprecated parameters for the initialization function; pass a single object instead")), r === void 0 && (r = new URL("pymol_web_bg.wasm", /* @__PURE__ */ ((c) => c)(import.meta.url)));
  const e = E();
  (typeof r == "string" || typeof Request == "function" && r instanceof Request || typeof URL == "function" && r instanceof URL) && (r = fetch(r));
  const { instance: t, module: _ } = await ye(await r, e);
  return $(t);
}
export {
  T as WebViewer,
  ve as default,
  he as init,
  Se as initSync
};
